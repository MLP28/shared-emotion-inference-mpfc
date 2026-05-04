# Socially shared emotions shape the activity of the medial prefrontal cortex during inference of others' emotional states

This repository contains the MATLAB analysis code associated with the paper:

> **Socially shared emotions shape the activity of the medial prefrontal cortex during inference of others' emotional states**  
> Le Petit, Gagnepain, ... Laisney, 2026 Cell Reports

---

## Overview

The code implements three fMRI analysis pipelines using SPM12:

- **Univariate analysis** — preprocessing, first-level GLM, and second-level group statistics
- **RSA (Representational Similarity Analysis)** — searchlight RSA using the [Viviani et al. toolbox](https://github.com/viviani-rsa/rsa-rsm)
- **PPI (Psychophysiological Interaction)** — seed-based functional connectivity analysis

---

## Repository Structure

```
matlab_codes_for_publication/
│
├── data/
│   ├── COST/                        # Raw MRI data (group 1)
│   ├── COSSA/                       # Raw MRI data (group 2)
│   ├── seed_masks/                  # ROI masks used as seed in PPI analysis
│   └── behav/
│       ├── trialinfo_MRIstudy/      # Subject-level onset files
│       ├── infRT_prestudy/          # Reaction times from pre-study 2
│       ├── emointensity_prestudy/   # Emotion ratings from pre-study 1
│       └── infotask/                # Task information (infoemo.xlsx)

│
├── dependency/
│   └── rsa-rsm-master/              # Viviani RSA toolbox (external)
│
├── fun/                             # Shared utility functions (get_files.m, etc.)
│
├── fun_univariate/                  
│   ├── cos_run_preproc_1lev.m       # Univariate entry function for preprocessing and first level
│   ├── cos_run_secondlevel.m        # Univariate entry function for second level
│   ├── fun/                         # univariate helper functions
│   └── error/                       # log in case of errors during univariate analysis
│
├── fun_RSA/                         # RSA-specific functions
│   ├── cos_rsa_run_firstlevel.m     # RSA entry function for first level
│   ├── cos_rsa_run_secondlevel.m    # RSA entry function for second level
│   ├── fun/                         # RSA helper functions
│   └── error/
│
├── fun_PPI/                         # PPI-specific functions
│   ├── cos_run_ppi.m                # PPI entry function for first level
│   ├── cos_run_ppi_secondlevel.m    # PPI entry function for second level 
│   ├── fun/                         # PPI helper functions
│   └── error/
│
├── results_univariate/                
│   ├── COST/                       # MRI data after univariate analyses (group 1)
│   ├── COSSA/                      # MRI data after univariate analyses (group 2)
│   └── SecondLevel/                # MRI data after second level univariate analyses
│
├── results_RSA/                
│   ├── COST/                       
│   ├── COSSA/ 
│   ├── rsm/                       # RSMs computed from pre-studies used in searchlight analysis
│   └── SecondLevel/   
│
├── results_univariate/                
│   ├── COST/                       
│   ├── COSSA/   
│   └── SecondLevel/  

```

---

## Dependencies

- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) — neuroimaging analysis toolbox
- [Viviani RSA toolbox](https://github.com/viviani-rsa/rsa-rsm) — searchlight RSA (Viviani et al., 2021)
- MATLAB Statistics and Machine Learning Toolbox (`fitdist`)

> **Important:** SPM orthogonalization must be manually disabled before running the univariate first-level analysis. Comment out line 116 in `spm_get_ons.m` and lines 242–244 in `spm_fMRI_design.m`. See [this reference](http://imaging.mrc-cbu.cam.ac.uk/imaging/ParametricModulations) for details.

---

## Usage

All analyses are controlled via **entry point scripts**. Open the relevant script, fill in the `USER SETTINGS` section at the top, and run.

### Univariate analysis

```matlab
% Preprocessing + first-level GLM (subject indices as input)
cos_run_preproc_1lev(1:N)

% Second-level group statistics
cos_run_secondlevel('one')    % one-sample t-test
cos_run_secondlevel('two')    % two-sample t-test
cos_run_secondlevel('both')   % both
```

### RSA

```matlab
% Build RDMs/RSMs from pre-study data (run once, before subject loop)
cos_rsa_run_firstlevel([], 'rdm')

% First-level pipeline per subject (reslice -> GLM -> searchlight -> postprocess)
cos_rsa_run_firstlevel(1:N, 'all')

% Individual steps can be run separately
cos_rsa_run_firstlevel(1:N, 'reslice')
cos_rsa_run_firstlevel(1:N, 'spm')
cos_rsa_run_firstlevel(1:N, 'searchlight')
cos_rsa_run_firstlevel(1:N, 'postprocess')

% Second-level group statistics
cos_rsa_run_secondlevel('one')
cos_rsa_run_secondlevel('two')
cos_rsa_run_secondlevel('both')
```

### PPI

```matlab
% First-level PPI maps per subject
cos_run_ppi(1:N)

% Second-level group statistics
cos_run_ppi_secondlevel('one')
cos_run_ppi_secondlevel('two')
cos_run_ppi_secondlevel('both')
```

---

## Configuration

Each entry point script contains a clearly delimited `USER SETTINGS` section at the top. Only this section needs to be edited before running. It includes:

- SPM and toolbox paths
- Root data and results directories
- Subject list and group names
- Analysis-specific parameters (contrasts, searchlight size, smoothing kernel, etc.)

Parameters that must be consistent across scripts (e.g., `glm_model`, `searchlight.size`, `searchlight.model_name`) are flagged with comments in the code.

---

## Output Structure

```
results_univariate/
    <group>/<subject>/
        run_*/          # preprocessed functional data
        anat/           # preprocessed anatomical data
        glm_emo/        # first-level SPM model
    SecondLevel/
        OneSample/<contrast>/
        TwoSample/<contrast>/

results_RSA/
    rsm/                # RSM files (model RDMs)
    <group>/<subject>/
        preprocessing/  # resliced functional images
        glm_rsanative/  # first-level RSA GLM
        searchlight/    # individual RSA maps
    SecondLevel/
        OneSample/<rsm_name>/
        TwoSample/<rsm_name>/

results_PPI/
    <subject>/<seed>/   # individual PPI maps
    SecondLevel/
        OneSample/<seed>/<contrast>/
        TwoSample/<seed>/<contrast>/
```

---

## Contact

Marine Le Petit (marine.lepetit@uliege.be)
Pierre Gagnepain (gagnepain@cyceron.fr)
Mickaël Laisney (mickael.laisney@ephe.psl.eu)
U1077 NIMH, Caen, France
