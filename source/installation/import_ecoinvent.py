import bw2data as bd
import bw2io as bi
import brightway2 as bw
from bw2io.migrations import create_core_migrations
from collections import defaultdict
import numpy as np

bd.projects.set_current(
    "ecoinvent311"
)  # changer le nom du projet pour importer ecoinvent dans un nouveau projet

# Bout de code adapté de la fonction bwio.ecoinvent.import_ecoinvent_release(), pour pouvoir importer la version téléchargée d'ecoinvent

ei_path = r"S:\370-Energie\370.25-BATTERIES\370.25.311-Formation_Brightway\ecoinvent 3.11_cutoff_ecoSpold02"  # changer le path si autre version d'ecoinvent voulue

# Import biosphere
eb = bi.importers.Ecospold2BiosphereImporter(
    name="ecoinvent-3.11-biosphere",
    filepath=ei_path + r"\MasterData\ElementaryExchanges.xml",
)
eb.apply_strategies()
eb.write_database()
bd.preferences["biosphere_database"] = "ecoinvent-3.11-biosphere"

# Import ecoinvent-3.11-cutoff. prend environ 15 min
ei = bi.SingleOutputEcospold2Importer(ei_path + r"\datasets", "ecoinvent-3.11-cutoff")
ei.apply_strategies()  # apply functions to normalize the flows, turn into a adequate format and match all the properties of the flows (name, unit, reference product, location...)
ei.statistics()  # to check if some flows have no match
ei.write_database()

# Import LCIA methods
biosphere_name = bd.config.biosphere

lcia_file = ei_path + r"\LCIA Implementation 3.11.xlsx"
sheet_names = bi.ecoinvent.get_excel_sheet_names(lcia_file)

if "units" in sheet_names:
    units_sheetname = "units"
elif "Indicators" in sheet_names:
    units_sheetname = "Indicators"
else:
    raise ValueError(f"Can't find worksheet for impact category units in {sheet_names}")

data = dict(bi.extractors.ExcelExtractor.extract(lcia_file))
units = bi.ecoinvent.header_dict(data[units_sheetname])
cfs = bi.ecoinvent.header_dict(data["CFs"])
CF_COLUMN_LABELS = {
    "3.4": "cf 3.4",
    "3.5": "cf 3.5",
    "3.6": "cf 3.6",
}
cf_col_label = CF_COLUMN_LABELS.get("3.11", "cf")
units_col_label = bi.ecoinvent.pick_a_unit_label_already(units[0])
units_mapping = {
    (row["method"], row["category"], row["indicator"]): row[units_col_label]
    for row in units
}
biosphere_mapping = {}
for flow in bd.Database(biosphere_name):
    biosphere_mapping[(flow["name"],) + tuple(flow["categories"])] = flow.key
    if flow["name"].startswith("[Deleted]"):
        biosphere_mapping[
            (flow["name"].replace("[Deleted]", ""),) + tuple(flow["categories"])
        ] = flow.ke

lcia_data_as_dict = defaultdict(list)
unmatched = set()
substituted = set()
for row in cfs:
    impact_category = (row["method"], row["category"], row["indicator"])
    if row[cf_col_label] is None:
        continue

    lcia_data_as_dict[impact_category].append(
        (
            biosphere_mapping[
                bi.ecoinvent.drop_unspecified(
                    row["name"], row["compartment"], row["subcompartment"]
                )
            ],
            float(row[cf_col_label]),
        )
    )
for key in lcia_data_as_dict:
    method = bd.Method(key)
    method.register(
        unit=units_mapping.get(key, "Unknown"),
        filepath=str(lcia_file),
        ecoinvent_version="3.11",
        database=biosphere_name,
    )
    method.write(lcia_data_as_dict[key])

# create core migrations
create_core_migrations()

# Calcul de vérification
# setup
CC = [
    m for m in bd.methods if "EF v3.1 no LT" in str(m) and "climate change" in str(m)
][0]
AP = [m for m in bd.methods if "EF v3.1 no LT" in str(m) and "acidification" in str(m)][
    0
]
PM = [
    m
    for m in bd.methods
    if "EF v3.1 no LT" in str(m) and "particulate matter" in str(m)
][0]
ADP = [
    m
    for m in bd.methods
    if "EF v3.1 no LT" in str(m) and "material resources" in str(m)
][0]
list_methods = [CC, AP, PM, ADP]

ei = bd.Database("ecoinvent-3.11-cutoff")
coalDE = [
    a
    for a in ei
    if "electricity production" in a["name"]
    and "coal" in a["name"]
    and a["location"] == "DE"
][0]

FU = [{coalDE: 1}]
bw.calculation_setups["setupname"] = {"inv": FU, "ia": list_methods}
mLCA = bw.MultiLCA("setupname")

# Check que les résultats de l'ACV sont bons
assert np.allclose(
    mLCA.results,
    np.array([[1.08098023e00, 3.21409401e-03, 1.04029559e-08, 2.46450105e-07]]),
)

print("Ecoinvent v3.11 imported successfully!")
