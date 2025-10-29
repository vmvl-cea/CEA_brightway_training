# I. Activity Browser

## Open project

### Open Activity-Browser

Use commands in the Miniforge prompt:
• Search Miniforge prompt and open it
• Run the following commands (replace conda with micromamba if you have it installed):
conda a c t i v a t e bw2
a c t i v i t y −browser
and wait few seconds. "bw2" is the environment’s name where AB is installed.

### Import data

A partial LCA modelisation has been created on AB for this TP. You will download source files and install them on AB.
Download source files:
S:\370-Energie\370.25-BATTERIES\370.25.311-Formation_Brightway\tuto_EV
Import project:
• In your AB interface, in Projects, select "ecoinvent311".
• Import Database:
– Import tuto_EV-ab.xslx by clicking on Import, then Import local data, then Local
Excel file, then Browse and search the file.
After clicking on Import Database, a new window DATABASE LINKING opens. Select ecoinvent-3.11-cutoff in the drop-down menu

## Compute and understand impacts

### Add a calculation setup

In LCA Setup, add a new calculation setup by clicking on New, and name it Lifecycle.

### Choose impacts categories

AB only computes impacts you choose. You then need to select the indicators you want to compute. You can choose in this list of EF 3.1 indicators:
• climate change | global warming potential (GWP100)
• material resources: metals/minerals | abiotic depletion potential (ADP): elements (ul- timate reserves)
• energy resources: non-renewable | abiotic depletion potential (ADP) : fossil fuels
• acidification | accumulated exceedance (AE)
• ecotoxicity: freshwater | comparative toxic unit for ecosystems (CTUe)
• human toxicity: carcinogenic | comparative toxic unit for human (CTUh)
• human toxicity : non-carcinogenic | comparative toxic unit for human (CTUh)
• eutrophication: freshwater | fraction of nutrients reaching freshwater end compartment (P)
• eutrophication: marine | fraction of nutrients reaching marine end compartment (N)
• eutrophication: terrestrial | accumulated exceedance (AE)
• ionising radiation: human health | human exposure efficiency relative to u235
• land use | soil quality index
• ozone depletion | ozone depletion potential (ODP)
• particulate matter formation | impact on human health
• photochemical oxidant formation: human health | tropospheric ozone concentration increase
• water use | user deprivation potential (deprivation-weighted water consumption)
For this TP, choose the Climate Change indicator + 1 or 2 other indicators by
adding them in your calculation setup in Impact categories by a click and drag from the Impact Categories tab

### Calculate and interpret the impacts

#### Lifecycle

• Open your excel imported database tuto_EV-ab, and search Lifecycle.
• Click on it to read the description of this activity, and drag it to Reference flows.
• Open your ecoinvent-3.11-cutoff db and search for the activity of an internal combustion engine (ICE) car called market for transport, passenger, car, petrol, medium size, EURO 4, drag it in Reference Flow to compare this activity with the Lifecycle of the EV vehicle one.
• Click on calculate.
• Check the Climate change results

Do you think the activity Lifecycle has a great/little/normal impact for an EV in France? and compared to the ICE car ?
Look at the Sankey Diagram for Lifecycle and identify an anomaly. You will correct it at the next part.
Graph Explorer:
Go to the tuto_EV-ab database, and right-click on Lifecycle. Select Open activity in Graph Explorer. Explore the graph.

#### Production – contribution

- Add a new calculation setup, and name it Production
- Drag and drop the activity Production from your tuto_EV-ab database to Reference flows
- Drag and drop the impact categories you want to calculate
- Click on Calculate
- In the “LCA Results” tab, go to FT Contributions
You can see that the battery accounts for a big part of the vehicle’s production impacts !

## Edit model

### Modify activities

Delete exchanges:
Uncheck the Read-only box in the database list, and go to the activity identified in section 2.3.1 to delete the anomaly
(right-click > delete Exchange).
Add new exchanges :
The electricity consumed during the lifetime of the battery in the vehicle needs to be added. Add it to Use, with a click and drag from ecoinvent-3.11-cutoff database to Technosphere flows, and change the amount:
Modifications
Modifications areare lockedlocked atat 22 levels.levels. DatabasesDatabases areare lockedlocked whenwhen ReadRead--onlyonly isis checked,checked, and and activitiesactivities areare lockedlocked whenwhen EditEdit activityactivity (in(in eacheach activity)activity) isis checked.checked.

### Add formulas

• In the Parameter tab in Database, create 3 "Database parameters", as in the picture below:
• Go to Use activity, and write the following formula:
• Go to Lifecycle activity, and write the following formulas:
• Calculate the impacts generated in Lifecycle calculation setup. Go to Overview (in LCA results tab) and export table in Excel or CSV format.
• Go back to Parameters and modify efficiency value to 100 (Wh/km). Calculate the new impacts and compare them with your export.
Is the efficiency of the vehicle a significant factor in LCA results? For all your indicators?

### Add scenarios

effivciency
efficiency

Use a pre-built scenario:
• In LCA setup tab, take the Lifecycle calculation setup and change Standard LCA to Scenario LCA.
• In Scenario, click on Add scenario and browse Scenarios-tuto_EV-ab.xlsx.
• Visualize the parameters values in Parameters tab > Scenarios tab.
• Calculate the impacts and interpret the results.
Build your own scenarios:
Modify the Excel file to create new scenarios:

1. Tesla model Y (<https://ev-database.org/car/3103/Tesla-Model-Y-RWD>)
2. Total_distance = 300 000 km
3. Whatever you want
Import the new Excel, calculate and interpret the results.
