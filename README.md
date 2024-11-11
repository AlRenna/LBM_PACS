# LBM_PACS


## SETUP
### Required packages
run 'sudo apt-get install "package"':
- ffmpeg
-

### Create the environment
```bash
python3 -m venv env
source env/bin/activate
pip install -r py_requirements.txt
```

### Create the build folder to execute the code
```bash
mkdir build
cd build
cmake ..
make
```