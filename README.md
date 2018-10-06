	_____________________________________________________________________________________________

				Sparsity-based super-resolution correlation microscopy (SPARCOM)
	_____________________________________________________________________________________________


** Contents

	1. Overview
	2. Requirements
	3. Installation and basic operation
	4. Copyright
	5. Warranty
	6. History
	7. Download
	8. Trademarks

** Publisher
	
	Oren Solomon				orensol@campus.technion.ac.il
		
	Department of Electrical Engineering
	Technion - Israel Institute of Technology
	Haifa, 32000, Israel

1. Overview

	We present a MATLAB code for sparsity-based super-resolution correlation microscopy (SPARCOM). 
	SPARCOM utilizes sparsity in the correlation domain and uncorrelated blinking statistics
	to perform super-resolution imaging of fluorescent microscopy sequences with high-emitter density. 
	The code provided here is self-contained and includes all necessary functions (except for the BM3D code, which can be downloaded at http://www.cs.tut.fi/~foi/GCF-BM3D/). 
	Currently, the implementation is single core, but can be extended to multi-core processing.
	
	The code includes an example configuration file and an example simulated movie, along with its by-products. The simulation dataset is freely available from http://bigwww.epfl.ch/smlm/#&panel1-1
	
	This code is for academic purposes only.
	
	If you are using this code, please cite: 
    1.	Solomon, Oren, et al. "Sparsity-based super-resolution microscopy from correlation information." Optics Express 26.14 (2018): 18238-18269.


2. Requirements

	• MATLAB R2016a or newer (previous versions require modifications).
	• At least 8GB RAM; 64GB or more is recommended.


3. Installation and basic operation
	
	To Install:
	-----------
	1. Unpack the archive and add the SPARCOM_V1 folder to the MATLAB path.
	2. AddPathList.m - adds all releveant directories to the path. Can be run once.
	
	To run    :
	-----------
	1. SPARCOM.m         - Main file. Just need to run it.
		1. VERBOSE         : 1 - Display massages to screen, 0 - do not display.
		2. InternalSaveFlag: 1 - save results to disk, 0 - do not save.
	2. SRFM_Config_1.txt - Configuration file. Parameters are changed from this file. Description for different parameters are inside the file.
	
	Main internal functions:
	------------------------
	'CoreFunctions\':
	1. InitialMovieProcessing.m - Loads the data and cuts it to the appropriate dimensions, as well as performing optional denoising.
	2. FeqPrep.m     	        - Preprocessing function which generates the Fourier domain implementation.
	3. PatchManipulator_V2.m    - Extracts low-resolution patches and combines high-resolution outputs.
	4. ReadConfigFile.m         - Reads the TXT configuration file.
	5. CorrMic.m				- Main processing script.
	6.CorrMic_Patches.m         - Main processing script, patch based analysis.
	
	'OptEngine\': Correlations solvers
	1. pFISTA_diag.m          - FISTA based l_1 recovery.
	2. pFISTA_iterative.m     - FISTA iterative l_1 reconstruction.
	3. pFISTA_diag_analysis.m - FISTA based analysis formulation recovery. Supports wavelet, DCT and fast Welch-Hadamard bases.
	4. pFISTA_diag_TV.m       - FISTA based TV recovery.
	5. pFISTA_BM3D.m      	  - FISTA based BM3D recovery.
	
	
4. Copyright

    Copyright © 2018 Oren Solomon, Department of Electrical Engineering, 
	Technion - Israel Institute of Technology, Haifa, 32000, Israel
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
	This code is for academic purposes only.
	
5. Warranty

	Any warranty is strictly refused and you cannot anticipate any financial or
	technical support in case of malfunction or damage; see the copyright notice.

	Feedback and comments are welcome. We will try to track reported problems and
	fix bugs.

	Bugs are encouraged to be reported to orensol@campus.technion.ac.il
	
6. History

  • August 12, 2018
	Version 1.0 released under GNU GPL version 3.


7. Download

	The code is available at http://webee.technion.ac.il/Sites/People/YoninaEldar/software.php
	and on GitHub: https://github.com/KrakenLeaf/SPARCOM 


8. Trademarks

	MATLAB is a registered trademark of The MathWorks. Other product or brand
	names are trademarks or registered trademarks of their respective holders.
	
	All third party software and code packages are disctributed under the GNU license as well. 
	The authors claim no responsibility for this software and code.
