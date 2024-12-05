### 0. Save environment recipes

#### 0-1. Save installed packages in YAML file

```
$ conda env export -n [environment_name] --no-build > [enviornment_name].yml
```

#### 0-2. Modify YAML file no to use anaconda channel when recreating it

- Remove `defaults` from channels in YAML file (.yml)

### 1. Remove miniconda (or anaconda)

#### 1-1. Remove init scripts in .bashrc

```
$ conda activavte
$ conda init --reverse --all	# Remove init script
$ conda clean --all		# Remove tarballs and caches
```

#### 1-2. Remove miniconda (or anaconda) directory

```
$ CONDA_BASE_ENVIRONMENT=$(conda info --base)

$ echo The next command will delete all files in ${CONDA_BASE_ENVIRONMENT}
# Warning, the rm command below is irreversible!
# check the output of the echo command above
# To make sure you are deleting the correct directory

$ rm -rf ${CONDA_BASE_ENVIRONMENT}
```

#### 1-3. Delete .conda/ and rename .condarc

```
$ rm -rf ~/.conda
$ mv ~/.condarc ~/.condarc.old
```

#### Log out, then in again

### 2. Install Miniforge3 (third-party counterpart of miniconda)

#### 2-1. Download install file from github and run

```
$ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
$ bash Miniforge3-$(uname)-$(uname -m).sh
# Follow the installation prompts
```

#### Log out, then in again

#### 2-2. (For ECBSU members) Add settings

```
$ conda config --append channels bioconda

# If you don't prefer auto activate of base environment
$ conda config --set auto-activate_base false
```

### 3. Recreate environments

#### Run the commands below for each enviroment
```
$ conda activate
$ conda env create -f [environment_name].yml
# Caution: When the package is available only in the anaconda channel, the installation will fail.
```


### References
- Miniforge (https://github.com/conda-forge/miniforge)
- Anaconda at OIST (https://groups.oist.jp/scs/anaconda-oist)