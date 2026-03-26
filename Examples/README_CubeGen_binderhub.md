# CubeGen Notebook — Data Cube Reconstruction Pipeline

This notebook provides an example of how to use 3DCubeGen within a notebook environment using BinderHub and provides a complete workflow to:

- Generate IFU data cubes  
- Produce emission-line maps  
- Create RGB composite images  
- Apply 2D deconvolution  
- Generate final science-quality products  

This workflow is designed to work with **LVM DRP outputs** and the **CubeGen / 3DCubeGen** package.

---

# 📚 Documentation

For full documentation of **3DCubeGen**, please visit:

https://hjibarram.github.io/3DCubeGen/

GitHub Repository:

https://github.com/hjibarram/3DCubeGen

---

# Requirements

```bash
git clone https://github.com/hjibarram/3DCubeGen.git
cd 3DCubeGen
pip install .
```

You also need:

- LVM DRP outputs

---

# Imports

```python
from CubeGen.tools.cubegen import map_ifu
from CubeGen.tools.mapgen import gen_map
from CubeGen.tools.images import get_jpg
from CubeGen.tools.deconvolve import deconvolve_2dfile
import CubeGen.tools.tools as tools
```

---

# Step 1 — Generate Data Cube

Generate IFU data cube using:

```python
map_ifu(...)
```

Output:

```
lvmCube-*.fits
```

---

# Step 2 — Generate Emission Line Maps

```python
gen_map(...)
```

Creates:

- OIII maps
- Hα maps
- SII maps

---

# Step 3 — Create RGB Composite Image

```python
get_jpg(...)
```

---

# Step 4 — Apply Deconvolution

```python
deconvolve_2dfile(...)
```

---

# Step 5 — Final RGB

```python
get_jpg(...)
```

---

# Workflow

Raw LVM DRP  
↓  
map_ifu  
↓  
gen_map  
↓  
get_jpg  
↓  
deconvolve_2dfile  
↓  
get_jpg  

---

# Output Files

- lvmCube-*.fits  
- lvmMap-*_OIII.fits  
- lvmMap-*_HI.fits  
- lvmMap-*_SII.fits  
- lvmMap-*_OHS.jpeg  
- lvmMap-*_decv.fits  
- lvmMap-*_OHS_decv.jpeg  

---

# Documentation

https://hjibarram.github.io/3DCubeGen/

---

# GitHub

https://github.com/hjibarram/3DCubeGen
