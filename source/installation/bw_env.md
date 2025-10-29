# Setting up a Python Environment for Brightway

## ⚙️ Managing Python Environments with Conda

**Conda** is a popular tool to create and manage these Python environments. It comes with the free [Miniforge](https://conda-forge.org/download/) installer (authorized for use at CEA), which installs the `conda` and `mamba` commands on your computer.

Note: [Anaconda](https://www.anaconda.com/download) and [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) **cannot** be used as the CEA does not pay a licence for them, cf. the [STIC notice](https://portail.intra.cea.fr/grenoble/sti/Pages/Vous%20aider/Poste%20de%20travail/Anaconda.aspx))

### 1. Installing Conda

Install `Miniforge` from CEA's software center.

Open a terminal (search for 'cmd' or 'command prompt' in Windows Start menu) and run the following line:

```cmd
conda init
```

This will activate the `conda` commands in your terminal.

Close the terminal.

### 2. Create a new environment with Brightway

Download the environment requirement file here: [bw_training.yml](filepath) **TODO** Add file path.

Open a new terminal and run the following line:

```cmd
conda create -n bw25 -f C:\<path_to_requirement_file>\.yml
```

This creates a new environment named `bw25` with Python version 3.11, Brightway 2.5, Activity-browser 3 and lca_algebraic_25.

Note : we recommand keeping this file so that you can easily delete and recreate a fresh and clean environment for brightway 2.5.

#### 2. Activate the environment

```cmd
conda activate brightway-env
```

Now you’re *inside* your environment — any library you install will stay in this environment

Note : we recommand keeping this file, so that you can easily delete and recreate a fresh and clean environment for brightway 2.5.

#### Recreating fresh and clean environment

In case something goes wrong in your environment, you can always delete it and recreate a new clean one. No worries, this will not delete any of your Brightway projects or any other data. This will just uninstall all the libraries and delete.

-> Delete an environment:
-> Recreate an environment

--
Footnote and sources:

* written with the help of ChatGPT
* <https://www.dataquest.io/blog/a-complete-guide-to-python-virtual-environments/>
* <https://medium.com/data-science/what-is-a-python-environment-for-beginners-7f06911cf01a>
* <https://realpython.com/python-virtual-environments-a-primer/>
