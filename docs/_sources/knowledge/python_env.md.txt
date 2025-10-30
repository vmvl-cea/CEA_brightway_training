# Python Environment for Brightway

To use Brightway on your computer, we install it in a specific **Python Envrionment**. This is a general Python good practice, allowing to keep a clean version of Brightway (and all the packages it requires) isolated from any other Python project you may have.

## 🐍 What is a Python Environment?

A **Python environment** is like a *workspace* where you can install and use specific Python tools and libraries without interfering with other projects on your computer.

Think of it as a **sandbox** that keeps everything isolated:

* Each environment have its **own version of Python** (e.g., 3.10 vs. 3.12), which can thus be different than the base one on your computer.
* Each environment have **its own set of packages** (e.g., Brightway, NumPy, Pandas, etc.).
* This prevents conflicts — for example, one project might need an older library version, while another needs the latest one.

This isolation is especially useful in scientific computing, where different tools often have strict package version requirements.

You will see later how to create and manage this environment.

On your computer, you can have as many python environment as you wish, in addition to the `base` one (the one used by default).

![Alt Text](../images/python_virtual_environments.png)

They are stored locally in your computer, typically here : `C:\Users\AB123456\AppData\Local\miniforge3\envs`. Inside this folder, there is the Python executable `python.exe` of the Python version of the environment. There is also the code of all the librairies and packages (in `env_name\lib\site-package`). This is why you can have a different version of Python and different versions of the libraries in each environment.

## Python Environment and Brightway Projects

**Python Environments** on your computers are just workspaces in which to run python files and python-based software. These files and folders are independant from the environments. Thus:

* Creating or deleting Python environments has no impact on other projects files and folders, e.g. Brightway files and folders
* Creating or deleting files and folders related to python-based software has no impact on Python environments.

Brightway projects files are each stored in project-folders on your computers, usually here: `C:\Users\AB123456\AppData\Local\pylca\Brightway3`
  
**TODO** Picture with bw files and folders compared to python env files and folders.
