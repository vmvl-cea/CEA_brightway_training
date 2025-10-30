# II. Brightway

## 1. Install Brightway in a dedicated environment

1. Download the [environment requirements file](./bw25_training.yml) and save it somewhere convenient.

2. Open a new terminal, go to the location where you saved the file (`cd` command) and run the following command:

    ```cmd
    conda create -n bw25 -f bw25_training.yml
    ```

This will create a new environment named `bw25` and install Brightway and other libraries, as specified in the YAML file:

* `Python=3.12`
* latest version of `Brightway25`
* latest version of `Activity-browser` compatible with Brightway 2.5
* latest version of `Lca-algebraic-25`
* a couple of other useful libraries

Note: You can name your environment differently by replacing `bw25` with the name of your choice.

## 2. Verify the environment

1. Activate the environment by running the conda `activate` command:

    ```cmd
    conda activate bw25
    ```

    You then see that the environment is active as its name should appear in parenthesis.

2. Verify the installation by running the `activity-browser` command, which should open Activity-Browser.

    ```cmd
    activity-browser
    ```

3. You can then deactivate your environment by running:

    ```cmd
    conda deactivate
    ```

## Recreating fresh and clean environment

In case something goes wrong in your environment, you can always delete it and recreate a clean new one. No worries, this will not delete any of your Brightway projects or any other data. This will just uninstall all the libraries and delete.

1. Delete your `bw25` environment by running:

    ```cmd
    conda env remove -n bw25
    ```

    Note: for the `remove` command to work, you must be outside of the environment.

2. Repeat the environment creation steps.
