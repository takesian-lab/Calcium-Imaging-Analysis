
INSTRUCTIONS FOR MAKING Z-CORRECTED BLOCKS


1) IMAGING WITH THE PIEZO
   Record continuous Z-series within a T-series BOT
        - 100000 reps (number of reps must be longer than Tosca file)
	- BOT checked
	- Max speed checked
	- Voltage Recording = current
	- Not 'In Use'
	- Set Piezo start and stop points to make a ~20um stack
        - Record 6 slices at 4um steps = 5fps
	- Leave shutter open
	- Fastest acquisition
	- One-directional
   Collect a Z-stack in the same position with greater precision and volume
        - 2um step size
        - Start and stop points don't matter as long as the volume includes the FOV


1.5) CONVERT BIDIRECTIONAL DATA TO ONE-DIRECTIONAL
   Suite2p expects multiplane data to be sequential. For now, I wrote a script to rename
   bidirectional tiffs to match what Suite2p expects.
   Run convert_bidirectional_tiffs


2) RUN DATA THROUGH SUITE2P
   Compiling with Suite2p is the same except for nplanes = 6 and fs = 5
   multiplane_parallel might make computation faster (I haven't tried yet)
   The data for each plane will be saved in separate folders (plane0, plane1... etc)
   You can compile multiple blocks together


3) PROCESS AND CROP Z-STACK FOR Z-CORRELATIONS
   The larger the Z-stack is, the longer it will take Suite2p to compute Z-correlations
   Use FIJI's slice remover to make the Z-stack smaller (i.e. reduce from 100 planes to 40)
   Save Z-stack as a tiff to upload into Suite2p

	Slice remover documentation: https://imagej.nih.gov/ij/plugins/slice-remover.html

	STEPS:
	1. Open z-stack in FIJI
	2. Load corresponding FOV in Suite2p (any plane)
	4. Find matching neurons in z-stack and note z plane of focal point
	5. Use slice remover to remove planes above and below focal point
	6. Reset brightness and save as tiff in block folder


4) COMPUTE ZCORR IN SUITE2P
   *Repeat this for every plane (i.e. plane0 stat.npy, then plane 1 stat.npy... etc)
   There is an option to load all planes at once (File, Load multiplane data) but I haven't tested this

	STEPS:
	1. Load stat.npy
	2. OPTIONAL: Define ROIs (can do either before or after Zcorr)
	3. View registered binary > load z-stack tiff > compute z position
	4. Track progress in terminal
	5. When finished processing, click save to mat. Check if zcorr is now saved in fall.mat.
	6. If not, reload stat.npy and try again (This is a bug with older versions of suite2p)
	

5) COMPILE BLOCKS
   Add blocks to info sheet
   *Direct analysis path (analysis_name) to suite2p folder instead of plane0 folder
   Run compile_blocks_from_info_V2
   *Message 'Found X suite2p planes' should appear to confirm that the block is being run as multiplane


6) VISUALIZE DATA
   Load block
   Run visualize_zcorr(block,plane)


7) MAKE ZCORRECTED BLOCKS
   Load block
   Run zcorrect_multiplane(block,1) -> plot figures only
       zcorrect_multiplane(block,0) -> make new block folder with tiffs from best planes only
                                       this will also make a new BOT with timestamps that correspond to the best planes
				       this will also copy voltage recording and XML files


8) RUN ZCORRECTED BLOCKS THROUGH SUITE2P, DEFINE ROIs, and COMPUTE ZCORR
   nplanes = 1, fs = 5
   You can compile multiple blocks together
   Define ROIs and compute Zcorr as in step 4


9) COMPILE ZCORRECTED BLOCKS
   Add blocks to info sheet. This should act like a normal block now.
   *Direct analysis path (analysis_name) to plane0 folder
   *block_name now starts with 'Zcorrected-'
   *IMPORTANT: Add '-Zcorrected' to the end of the stim_name (e.g. FMsweep-Zcorrected)
    otherwise it will save over your multiplane blocks because they have the same name


10) VISUALIZE DATA
   Load Zcorrected block
   Run visualize_zcorr(block)
   Admire lack of movement artifacts


11) OPTIONAL: SPLIT TIFFS INTO DIFFERENT FOLDERS BY PLANE
   If you want to split block tiffs by plane, you can use code split_BOTs_by_plane
   which will make folders plane0 - planeX and copy the corresponding tiffs into
   each folder and rename them so that you can make videos in FIJI or run them 
   through suite2p. This currently doesn't copy BOT, voltage recording or XML files
