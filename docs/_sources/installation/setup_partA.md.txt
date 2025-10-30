# I. Conda and VScode

In this section, we install Brightway 2.5 in a dedicated environment on your computer.

We will do so using the Python environment manager `conda`.

## 1. Install Conda

**Conda** is a popular tool to create and manage Python environments. It comes with the free [Miniforge](https://conda-forge.org/download/) distribution (authorized for use at CEA), which installs both `conda` and `mamba` commands in your terminals.

Note: [Anaconda](https://www.anaconda.com/download) and [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) **cannot** be used as the CEA does not pay a licence for them, cf. the [STIC notice](https://portail.intra.cea.fr/grenoble/sti/Pages/Vous%20aider/Poste%20de%20travail/Anaconda.aspx))

1. Install **Miniforge** from CEA's software center.

2. Open a terminal and activate **Conda** commands on your computer by running:

    ```cmd
    conda init
    ```

3. Verify that Conda is well installed. The following command should just print its version number:

    ```cmd
    conda -V
    ```

## 2. Install Visual Studio Code and extensions

Visual Studio Code (or *VS Code*) is a free, open-source code editor developed by Microsoft. It offers powerful features for code writing, debugging, Git integration, and a rich extension ecosystem to add functionnalities.

1. Install **VS code** from CEA's software center. Open it when the installation is finished

2. On the left sidebar, click the Extensions icon (looks like four squares).

3. Search for the official 'Python' extension from Microsoft and install it.

4. Repeat with the official Jupyter extension (also by Microsoft).
