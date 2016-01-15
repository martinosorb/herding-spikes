using System;
using System.Runtime.InteropServices;
using BW;
using System.IO;
using System.Collections;
using System.Windows.Forms;
using System.Diagnostics;
using System.Threading;
using System.Drawing;

//want an averaging over a longer interval for baseline estimate. make use of the asymmetric baseline which makes it easier to find an estimate
// iteratively take second lowest value over blocks of 4 frames

/*
Changes:
-Possibility to run the detection with reversed time axis, and/or to downsample the data by a factor of 2

-implemented additional filter for low frequencies (50Hz). This version 
averages blocks of 2 frames (except for downsampling case)
takes the second lowest value in blocks of 4 averages (done every 8 frames)
and again the second lowest value from 4 such blocks (total 32 frames, done every 16(<=12kHz) or 32 frames).
Then take a weighted sliding average of 3 (NFblocks, increase to 4 for >20kHz) of those estimates, 
where the most recent estimate increases its weight by 1/4 every 8 frames, and the weight of the oldest estimate is decreased by 1/4.
This procedure is done ahead of the detection time (TauFiltOffset),
and there is a parameter for the timescale to adapt a baseline estimate (TauFilt).

-converted matrices of surrounding electrodes into (seperate) jagged arrays, can use 'foreach' and seems to be a bit faster. Also 
made a list of locations that do not include locations at the boundary (can remove a check at runtime)

-Allows to keep a voltage mimimum for 3 frames (have to set dtEmx to 3 for that) 
to allow for larger temporal differences of the same spike on neighboring electrodes

-Compute weighted averages and test for threshold crossings only every dtTmx (i.e. 3) frames to make the algorithm faster.
when (Threshold-2) is crossed, every frame is analyzed. Further, the history of filtered voltages is kept for one frame 
and, when (Threshold-2) is crossed weighted averages are computed for the preceeding frame as well. 
(looses < 1% of detections, possibly with short spike shapes and mostly located on electrodes, therefore not necessarily real spikes)
Should be fine, since voltages are averaged over 3 frames and minima are kept for (dtEMx-1) frame 
(but the gap is effectively only (dtEMx-1) frames, and that for a lower threshold).

-Allows to define the cut-out region (excluding the center frame).

-Scaling amplitude and threshold initially rather than during the detection. Put decimal numbers in the gui 
(need to divide previous valued by 2; i.e. 5 instead of 10), so now they represent multiples of the variability estimate v.

- removed determination of correlations (have to put that back some time)

- have a parameter for the smoothing window (2 frames+Sampling/5000) for the detection (adjustable). 
This is one frame longer for the localization.

- Moved part of the localization here to save space in the output files and computation time in the postprocessing. Only output 4-5 (filtered) raw data traces (instead of 9-12).
Those will be used for determination of shapes (the others would have low weights anyway).
*/


namespace SpkDslowFilter {
	class Detection {
		//channels and their neighbors
		int NChannels;// number of channels; is set when reading the data
		int[][] ChInd4a;
		int[][] ChInd4b;
		int[][] ChInd4c;
		//int[][] ChInd4d;
		int[] ChInd4List;
		int[] ChInd5;
		int[][] ChInd5a;
		int[][] ChInd5b;
		int[][] ChInd5c;
		//int[][] ChInd5d;
		int[] ChInd5List;
		//Sampling rate and subsampling
		int Sampling;//=7000; set at runtime
		decimal FrameDt;//interval between frames. for gui.
		int sfi;// integer sampling rate
		bool HF;//frequency above 12kHz? Used for slow filter.
		bool Subsample;//whether only every second frame shall be used (used for slow filter).
		int[] dfTI;//to pass parameters (increment and offsets) from the gui to the main function.
		int df;//increment for detection (used for backward detection subsampling), defined in initial estimation
		int dfAbs;//absolute value (of df)
		int dfSign;//sign (of df)
		long t0x;//where the first cutout starts
		int tx;//where to start the iteration initially (relative to cutout)
		int ti;//initial point to estimate (relative to cutout)
		int tf;//final point to estimate (relative to cutout)
		int tm;//where to start Aglobal estimation
		int tms;
		int tfiA;
		int ty;//where to end Aglobal estimation
		//slow filter
		int[,] Avgs1;
		int[,] Avgs2;
		int[] Avgs3;
		int[] vSmth;
		int tSmth;
		int tSmth1;
		int tdSmth;
		int tdSmth1;
		const int TauFilt1 = 4;//i.e. 4*2=8 frames
		int TauFiltOffset;// = 48;//16 for smooth update (Taufilt=32) and 32 for overlapping, 64 for nonoverlapping. set at runtime.
		int NFblocks = 3; //average over two, weighted; want ~5 ms take Sampling/5000; set at runtime if >20kHz
		const int TauFilt=6*TauFilt1;//update timescale of slow offset voltage QmPreD in frames (here using 2*TauFilt1*TauFilt2)
		int FiltNorm;// = 18; set at runtime
		const int dtTMx = 3;
		//global fluctuations
		int[] Aglobal;
		int[] Aglobaldiff;
		int[] AglobalSmth;
		int[] AglobalSdiff;
		//int Aglobaldiffold=0;
		int[] Slice;
		//correlation with global fluctuations
		long SqIglobal;//sum of squared global increments
		long[] SqIv;//sum of squared channel increments
		long[,] SIprod;//sum of product of global and channel voltage increments
		int[] Vbias;
		int[] FVbias;
		long[] FVsbias;
		long[] Vsqbias;
		bool ACF;//measure autocorrelation? (otherwise faster)
		//Variables for variance and mean
		int[] Qd;//noise amplitude
		int[] Qm;//median
		int[] QmPre;//filtered raw
		int[] QmPreD;//offset after slow timescale correction
		int[,] Qdiff;//signal amplitude
		int[,] Qmax;
		int[] QmaxE;
		//Variables for the spike detection
		int[] Sl4;//counter for spike length
		int[] Z4next;//counter, skip frames when far from detection threshold
		int[] Z5next;
		bool[] Sl4x;//tag for removed spikes
		bool[] AHP4;//tag for repolarizing current
		int[] Amp4;//buffers spike amplitude
		int NSpikes4=0;//to count spikes
		//Variables for the spike detection
		int[] Sl5;//counter for spike length
		bool[] Sl5x;//tag for removed spikes
		bool[] AHP5;//tag for repolarizing current
		int[] Amp5;//buffers spike amplitude
		int NSpikes5=0;//to count spikes
		//cutouts
		int CutPre;
		int CutPost;
		int CutOffset;// = 13; set at runtime
		int tCut;//cutout interval=20; set at runtime
		int tCutLong;//=27; set at runtime
		int[] tQm;
		int[] tQmLong;
		int[] tQmA;
		int[] tQmX;//cutout indices
		int[] tQmXLong;//cutout indices
		int[] tShape;
		int tQmi;//initial cutout
		int tQmf;//final cutout
		int tQm0;//
		int tQmm;//median
		int Lspike;
		int Lsw;
		int Li;
		int Lf;
		int Lm;
		int[] Lw;
		int Lmx;
		//Parameters for variance and mean updates
		const int Tau_m0 = 4;//timescale for updating Qm (increment is Qd/Tau_m)
		const int Qdmin=300;//set minimum value of Qd
		//Parameters for spike detection
		int threshold;// = 6*AmpScale;//threshold to detect spikes >6 is likely to be real spikes; set at runtime
		int AHPthr;// = 0;//signal should go below that threshold within Slmax-Slmin frames; set at runtime
		int Slmax;// = 8;//dead time in frames after peak, used for further testing; set at runtime
		int Sln0;//=2; set at runtime
		//Parameters for reading data
		const int tInc = 192;//increment for reading data
		const int Ascale0 = -192;//factor to multiply to raw traces to increase resolution; definition of ADC counts had been changed!
		int Ascale;
		int AscaleV;
		int AscaleG;
		const int AmpScale=100;//want to use integer amplitudes
		//Parameters for recalibration events and artefact handling
		const int artT =10;//to use after artefacts; to update Qm for 10 frames
		int[] A;//control parameter for amplifier effects
		//instead of copying data around, want lookup variables that loop over one dimension of a matrix
		//2-loop
		int dt=0;
		const int dtMx = 2;
		int dtPre=1;
		//(dtEMx-1)-loop --> still exists/necessary?
		int dtE=0;
		const int dtEMx = 2;
		int dtEx=dtEMx-1;
		//Files to save the spikes etc.
		static FileStream fs;//for spikes
		StreamWriter w;
		static FileStream fsShapes;//for raw data
		StreamWriter wShapes;
		static FileStream fsX;//for spikes
		StreamWriter wX;
		static FileStream fsShapesX;//for raw data
		StreamWriter wShapesX;
		static FileStream fsInfo;//for other stuff
		StreamWriter wInfo;
		static FileStream fsMean;//for avg. Voltage
		StreamWriter wMean;
		int recalibTrigger=1;
		int Acal=3000;//for recalibration events
		// gui elements
		Form logwindow;
		TextBox logbox;
		/*
		//for testing only
		//inserted spikes
		int[,] Iweights = new int[16, 9] {{1000, 190, 250, 250, 250, 250, 190, 190, 190}, 
			{1000, 208, 250, 300, 250, 214, 174, 208, 174}, 
			{750, 225, 248, 375, 248, 187, 159, 225, 159}, 
			{500, 240, 240, 500, 240, 166, 146, 240, 146}, 
			{1000, 208, 300, 250, 214, 250, 208, 174, 174}, 
			{1000, 232, 300, 300, 214, 214, 187, 187, 161}, 
			{750, 258, 297, 375, 213, 187, 168, 198, 149}, 
			{500, 282, 282, 500, 208, 166, 153, 208, 138}, 
			{750, 225, 375, 248, 187, 248, 225, 159, 159}, 
			{750, 258, 375, 297, 187, 213, 198, 168, 149}, 
			{679, 297, 370, 370, 187, 187, 177, 177, 140}, 
			{486, 339, 339, 486, 183, 166, 159, 183, 131}, 
			{500, 240, 500, 240, 166, 240, 240, 146, 146}, 
			{500, 282, 500, 282, 166, 208, 208, 153, 138}, 
			{486, 339, 486, 339, 166, 183, 183, 159, 131}, 
			{414, 414, 414, 414, 163, 163, 163, 163, 123}};
		int[] Inn = new int[9] {0, 65, 64, 1, -64, -1, 63, -63, -65};
		int[] Ann = new int[15] {12, 11, 13, 10, 14, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19};
		int[] Iseq = new int[256] {65, 833, 1601, 2369, 3137, 3905, 577, 1345, 2145, 2913, 3681, 353, 1121, 1889, 2657, 3425, 69, 837, 1605, 2373, 3141, 3909, 581, 1349, 2149, 2917, 3685, 357, 1125, 1893, 2661, 3429, 73, 841, 1609, 2377, 3145, 3913, 585, 1353, 2153, 2921, 3689, 361, 1129, 1897, 2665, 3433, 77, 845, 1613, 2381, 3149, 3917, 589, 1357, 2157, 2925, 3693, 365, 1133, 1901, 2669, 3437, 81, 849, 1617, 2385, 3153, 3921, 593, 1361, 2161, 2929, 3697, 369, 1137, 1905, 2673, 3441, 85, 853, 1621, 2389, 3157, 3925, 597, 1365, 2165, 2933, 3701, 373, 1141, 1909, 2677, 3445, 89, 857, 1625, 2393, 3161, 3929, 601, 1369, 2169, 2937, 3705, 377, 1145, 1913, 2681, 3449, 93, 861, 1629, 2397, 3165, 3933, 605, 1373, 2173, 2941, 3709, 381, 1149, 1917, 2685, 3453, 97, 865, 1633, 2401, 3169, 3937, 609, 1377, 2113, 2881, 3649, 321, 1089, 1857, 2625, 3393, 101, 869, 1637, 2405, 3173, 3941, 613, 1381, 2117, 2885, 3653, 325, 1093, 1861, 2629, 3397, 105, 873, 1641, 2409, 3177, 3945, 617, 1385, 2121, 2889, 3657, 329, 1097, 1865, 2633, 3401, 109, 877, 1645, 2413, 3181, 3949, 621, 1389, 2125, 2893, 3661, 333, 1101, 1869, 2637, 3405, 113, 881, 1649, 2417, 3185, 3953, 625, 1393, 2129, 2897, 3665, 337, 1105, 1873, 2641, 3409, 117, 885, 1653, 2421, 3189, 3957, 629, 1397, 2133, 2901, 3669, 341, 1109, 1877, 2645, 3413, 121, 889, 1657, 2425, 3193, 3961, 633, 1401, 2137, 2905, 3673, 345, 1113, 1881, 2649, 3417, 125, 893, 1661, 2429, 3197, 3965, 637, 1405, 2141, 2909, 3677, 349, 1117, 1885, 2653, 3421};
		int[,] Template = new int[16,24] {{-406, 24694, 57565, 76477, 55086, 17420, -9309, -21894, -25392, -24455, -21813, -18796, -15964, -13508, -11454, -9764, -8385, -7260, -6341, -5587, -4966, -4450, -4021, -3660},
			{-146, 26462, 58798, 76391, 52769, 15302, -10468, -22316, -25420, -24315, -21622, -18606, -15796, -13365, -11335, -9668, -8306, -7196, -6288, -5544, -4929, -4420, -3996, -3639},
			{529, 28533, 60276, 76153, 50443, 13249, -11577, -22720, -25449, -24185, -21444, -18430, -15639, -13233, -11226, -9578, -8233, -7137, -6240, -5504, -4897, -4393, -3973, -3620},
			{1554, 30833, 61925, 75757, 48104, 11256, -12641, -23106, -25480, -24066, -21279, -18265, -15493, -13110, -11125, -9496, -8167, -7083, -6196, -5468, -4868, -4370, -3954, -3604},
			{2867, 33294, 63676, 75191, 45750, 9321, -13657, -23470, -25505, -23948, -21120, -18108, -15354, -12993, -11029, -9418, -8104, -7032, -6155, -5435, -4841, -4347, -3935, -3589},
			{4416, 35850, 65466, 74442, 43377, 7440, -14624, -23805, -25518, -23828, -20961, -17952, -15217, -12878, -10935, -9342, -8042, -6982, -6115, -5403, -4814, -4326, -3918, -3575},
			{6150, 38442, 67235, 73499, 40984, 5613, -15542, -24112, -25516, -23700, -20799, -17796, -15079, -12763, -10840, -9265, -7980, -6932, -6075, -5370, -4788, -4304, -3900, -3560},
			{8023, 41017, 68929, 72358, 38573, 3842, -16408, -24386, -25496, -23563, -20631, -17635, -14939, -12645, -10744, -9187, -7917, -6881, -6033, -5336, -4760, -4281, -3881, -3545},
			{9990, 43528, 70507, 71020, 36146, 2128, -17219, -24621, -25451, -23410, -20453, -17467, -14793, -12523, -10644, -9106, -7851, -6827, -5989, -5300, -4731, -4257, -3861, -3528},
			{12009, 45933, 71930, 69492, 33715, 475, -17974, -24819, -25382, -23241, -20263, -17290, -14640, -12396, -10539, -9020, -7782, -6771, -5943, -5262, -4699, -4230, -3839, -3509},
			{14042, 48201, 73174, 67788, 31284, -1116, -18674, -24979, -25289, -23056, -20063, -17105, -14481, -12263, -10430, -8931, -7709, -6711, -5893, -5221, -4665, -4202, -3815, -3489},
			{16060, 50301, 74222, 65927, 28873, -2641, -19319, -25106, -25176, -22859, -19855, -16915, -14318, -12127, -10318, -8839, -7633, -6649, -5842, -5179, -4629, -4172, -3789, -3467},
			{18032, 52219, 75066, 63930, 26488, -4103, -19916, -25202, -25046, -22652, -19641, -16721, -14152, -11989, -10204, -8746, -7557, -6586, -5790, -5135, -4593, -4141, -3763, -3444},
			{19936, 53946, 75703, 61822, 24141, -5497, -20465, -25272, -24903, -22439, -19423, -16525, -13985, -11850, -10090, -8652, -7479, -6522, -5737, -5091, -4555, -4110, -3736, -3421},
			{21756, 55477, 76142, 59629, 21842, -6827, -20973, -25322, -24752, -22224, -19206, -16331, -13819, -11712, -9976, -8559, -7403, -6459, -5684, -5047, -4518, -4078, -3709, -3398},
			{23479, 56814, 76391, 57372, 19598, -8096, -21445, -25358, -24599, -22012, -18994, -16141, -13658, -11578, -9866, -8469, -7329, -6397, -5634, -5004, -4483, -4048, -3684, -3376}};
		*/
		
		public int[] SetInitialParams (long nFrames, double nSec, int sf, double sfd, int NCh, int[] Indices)
		{
			dfTI= new int[3];
			Aglobal=new int[tInc];
			Aglobaldiff=new int[tInc];
			AglobalSdiff=new int[tInc];
			for (int i=0; i<tInc; i++) {
				Aglobal [i] = 4094;
				Aglobaldiff [i] = 0;
				AglobalSdiff [i] = 0;
			}
			//for (int i=0; i<2*dtMx-1;i++){
			//	qInd[i]=i%dtMx;
			//}
			//int[][] foo = new int[2][];
			//foo[0]=new int[] {10,11};
			//foo[1]=new int[] {3,4,5};
			//foreach (int i in foo[1]) {
			//	Console.WriteLine (i);
			//}


			NChannels = NCh;
			//NFrames = nFrames;
			sfi = sf / 1670;
			FrameDt = (decimal)(1000.0 / sfd);
			int[] SInd = new int[4096];
			int[] SInd4 = new int[4096];
			int[] SInd5 = new int[4096];
			ChInd4a = new int[NCh][];
			ChInd4b = new int[NCh][];
			ChInd4c = new int[NCh][];
			//ChInd4d = new int[NCh][];
			ChInd5 = new int[NCh];
			ChInd5a = new int[NCh][];
			ChInd5b = new int[NCh][];
			ChInd5c = new int[NCh][];
			//ChInd5d = new int[NCh][];
			for (int i=0; i<NCh; i++) {//fillvalues
				ChInd4a[i]= new int[4];
				ChInd4b[i]= new int[8];
				ChInd4c[i]= new int[4];
				//ChInd4d[i]= new int[4];
				ChInd5a[i]= new int[4];
				ChInd5b[i]= new int[4];
				ChInd5c[i]= new int[4];
				//ChInd5d[i]= new int[12];
				SInd[i]=-1;
				SInd4[i]=0;
				SInd5[i]=0;
				for (int j=0; j<8; j++) {
					ChInd4b [i][j] = -1;
				}
				for (int j=0; j<4; j++) {
					ChInd4a [i][j] = -1;
					ChInd4c [i][j] = -1;
					ChInd5a [i][j] = -1;
					ChInd5b [i][j] = -1;
					ChInd5c [i][j] = -1;
				}
				ChInd5 [i] = -1;
			}
			for (int i=0; i<NCh; i++) {//find active channels and number of neighbors
				SInd [Indices [i]] = i;
				SInd4 [Indices [i]] += 1;
				SInd5 [Indices [i]] += 1;
				if ((Indices [i] % 64) >0) {//channel to the left
					SInd4 [Indices [i] - 1] += 1;
					SInd5 [Indices [i] - 1] += 1;
					if ((Indices [i] / 64) >0) {//up and left
						SInd4 [Indices [i] - 65] += 1;
					}
				}
				if ((Indices [i] % 64) <63) {//right
					SInd5 [Indices [i] + 1] += 1;
				}
				if ((Indices [i] / 64) >0) {//up
					SInd4 [Indices [i] - 64] += 1;
					SInd5 [Indices [i] - 64] += 1;
				}
				if ((Indices [i] / 64) <63) {//down
					SInd5 [Indices [i] + 64] += 1;
				}
			}//channels with (all) neighbors have SInd4=4 or SInd5=5
			int ChInd4aN = 0;
			int ChInd5N = 0;
			for (int i=0; i<4096; i++) {
				if (SInd4 [i] == 4) {
					ChInd4a [SInd [i]][0] = SInd [i];
					ChInd4aN++;
					ChInd4a [SInd [i]][1] = SInd [i + 1];
					ChInd4a [SInd [i]][2] = SInd [i + 65];
					ChInd4a [SInd [i]][3] = SInd [i + 64];
					if (SInd5 [i] == 5) {
						ChInd4c [SInd [i]][0] = SInd [i];
					}
					if (SInd5 [(i + 1)%4096] == 5) {
						ChInd4c [SInd [i]][1] = SInd [(i + 1) % 4096];
					}
					if (SInd5 [(i + 64) % 4096] == 5) {
						ChInd4c [SInd [i]][3] = SInd [(i +64) % 4096];
					}
					if (SInd5 [(i + 65) % 4096] == 5) {
						ChInd4c [SInd [i]][2] = SInd [(i + 65) % 4096];
					}
					if (SInd4 [(i + 4032) % 4096] == 4) {
						ChInd4b [SInd [i]][0] = SInd [(i + 4032) % 4096];
						ChInd4b [SInd [i]][1] = SInd [(i + 4033) % 4096];
					}
					if (SInd4 [(i + 1) % 4096] == 4) {
						ChInd4b [SInd [i]][2] = SInd [(i + 2) % 4096];
						ChInd4b [SInd [i]][3] = SInd [(i + 66) % 4096];
					}
					if (SInd4 [(i + 64) % 4096] == 4) {
						ChInd4b [SInd [i]][4] = SInd [(i + 129) % 4096];
						ChInd4b [SInd [i]][5] = SInd [(i + 128) % 4096];
					}
					if (SInd4 [(i + 4095) % 4096] == 4) {
						ChInd4b [SInd [i]][6] = SInd [(i + 63) % 4096];
						ChInd4b [SInd [i]][7] = SInd [(i + 4095) % 4096];
					}
				}
				if (SInd5 [i] == 5) {
					ChInd5 [SInd [i]] = SInd [i];
					ChInd5N++;
					ChInd5a [SInd [i]][0] = SInd [i - 64];
					ChInd5a [SInd [i]][1] = SInd [i + 1];
					ChInd5a [SInd [i]][2] = SInd [i + 64];
					ChInd5a [SInd [i]][3] = SInd [i - 1];
					if (SInd4 [i] == 4) {
						ChInd5c [SInd [i]][2] = SInd [i];
						ChInd5b [SInd [i]][2] = SInd [(i + 65) % 4096];
					}
					if (SInd4 [(i + 4095)%4096] == 4) {
						ChInd5c [SInd [i]][3] = SInd [(i + 4095) % 4096];
						ChInd5b [SInd [i]][3] = SInd [(i + 63) % 4096];
					}
					if (SInd4 [(i + 4032) % 4096] == 4) {
						ChInd5c [SInd [i]][1] = SInd [(i + 4032) % 4096];
						ChInd5b [SInd [i]][1] = SInd [(i + 4033) % 4096];
					}
					if (SInd4 [(i + 4031)%4096] == 4) {
						ChInd5c [SInd [i]][0] = SInd [(i + 4031)%4096];
						ChInd5b [SInd [i]][0] = SInd [(i + 4031)%4096];
					}
				}
			}
			//ChInd5[ChInd4[i,ii],ii+5]

			//count -1 in ChInd4a[i][0]
			ChInd4List= new int[ChInd4aN];
			int iiii = 0;
			for (int i=0; i<NCh; i++) {
				if (ChInd4a [i] [0] != -1) {
					ChInd4List [iiii] = ChInd4a [i] [0];
					iiii++;
				}
			}
			ChInd5List= new int[ChInd5N];
			iiii = 0;
			for (int i=0; i<NCh; i++) {
				if (ChInd5 [i] != -1) {
					ChInd5List [iiii] = ChInd5 [i];
					iiii++;
				}
			}
			Qd = new int[NChannels];//noise amplitude
			Qm = new int[NChannels];//median
			QmPre = new int[NChannels];
			QmPreD = new int[NChannels];
			Qdiff = new int[NChannels,dtMx];
			Qmax = new int[NChannels,2];
			QmaxE = new int[NChannels];
			SqIv = new long[NChannels];//sum of squared channel increments
			SIprod = new long[NChannels,13];//sum of product of global and channel voltage increments
			//SIp = new int[NChannels];
			Vbias = new int[NChannels];
			FVbias = new int[NChannels];
			FVsbias = new long[NChannels];
			Vsqbias = new long[NChannels];
			A = new int[NChannels];//control parameter for amplifier effects
			Sl4 = new int[NChannels];//counter for spike length
			Sl4x = new bool[NChannels];
			Z4next = new int[NChannels];
			Z5next = new int[NChannels];
			AHP4 = new bool[NChannels];//counter for repolarizing current
			Amp4 = new int[NChannels];//buffers spike amplitude
			//Qd5 = new int[NChannels];//noise amplitude
			Sl5 = new int[NChannels];//counter for spike length
			Sl5x = new bool[NChannels];
			AHP5 = new bool[NChannels];//counter for repolarizing current
			Amp5 = new int[NChannels];//buffers spike amplitude
			Slice = new int[NChannels];
			Avgs1= new int[NChannels,3];
			//sortAvg= new int[NChannels,2];
			//Avgs1b= new int[NChannels][];
			Avgs3= new int[NChannels];
			// ask for some paramters
			Form form1 = new Form ();
			Button button1 = new Button ();
			button1.Text = "Go!";
			button1.Location = new Point (10, 10);
			button1.DialogResult = DialogResult.OK;
			form1.Text = "Detection Parameters";
			form1.FormBorderStyle = FormBorderStyle.FixedDialog;
			form1.AcceptButton = button1;
			form1.StartPosition = FormStartPosition.CenterScreen;
			form1.Controls.Add (button1);
			form1.MinimumSize= new Size(200,540);
			NumericUpDown numericUpDown1 = new NumericUpDown ();
			numericUpDown1.Parent = form1;
			numericUpDown1.Maximum = 10m;
			numericUpDown1.Minimum = 4m;
			numericUpDown1.Increment = 0.5m;
			numericUpDown1.Value = 6m;
			numericUpDown1.DecimalPlaces=1;
			numericUpDown1.Location = new Point (10, 70);
			Label l1 = new Label ();
			l1.Text = "Detection threshold";
			l1.AutoSize = true;
			l1.Location = new Point (numericUpDown1.Left, numericUpDown1.Bottom);
			l1.Size = numericUpDown1.Size;
			l1.Parent = form1;
			NumericUpDown numericUpDown3 = new NumericUpDown ();
			numericUpDown3.Parent = form1;
			numericUpDown3.Maximum = 2;
			numericUpDown3.Minimum = -2;
			numericUpDown3.Value = 1;
			numericUpDown3.Location = new Point (10, 220);
			Label l3 = new Label ();
			l3.Text = "Increment";
			l3.AutoSize = true;
			l3.Location = new Point (numericUpDown3.Left, numericUpDown3.Bottom);
			l3.Size = numericUpDown3.Size;
			l3.Parent = form1;
			NumericUpDown numericUpDown4 = new NumericUpDown ();
			numericUpDown4.Parent = form1;
			numericUpDown4.Maximum = 2.0m;
			numericUpDown4.Minimum = 0.0m;
			numericUpDown4.Increment = FrameDt;
			numericUpDown4.Value = 2*sfi*FrameDt;
			numericUpDown4.DecimalPlaces=2;
			numericUpDown4.Location = new Point (10, 270);
			Label l4 = new Label ();
			l4.Text = "CutoutPrePeak";
			l4.AutoSize = true;
			l4.Location = new Point (numericUpDown4.Left, numericUpDown4.Bottom);
			l4.Size = numericUpDown4.Size;
			l4.Parent = form1;
			NumericUpDown numericUpDown5 = new NumericUpDown ();
			numericUpDown5.Parent = form1;
			numericUpDown5.Maximum = 3.0m;
			numericUpDown5.Minimum = 0.0m;
			numericUpDown5.Increment = FrameDt;
			numericUpDown5.Value = 3*sfi*FrameDt;
			numericUpDown5.DecimalPlaces=2;
			numericUpDown5.Location = new Point (10, 320);
			Label l5 = new Label ();
			l5.Text = "CutoutPostPeak";
			l5.AutoSize = true;
			l5.Location = new Point (numericUpDown5.Left, numericUpDown5.Bottom);
			l5.Size = numericUpDown5.Size;
			l5.Parent = form1;
			NumericUpDown numericUpDown7 = new NumericUpDown ();
			numericUpDown7.Parent = form1;
			numericUpDown7.Maximum = 10*FrameDt;
			numericUpDown7.Minimum = 3*FrameDt;
			numericUpDown7.Increment = FrameDt;
			numericUpDown7.Value = (sf/5000+2)*FrameDt;
			numericUpDown7.DecimalPlaces=2;
			numericUpDown7.Location = new Point (10, 370);
			Label l7 = new Label ();
			l7.Text = "Smoothing kernel";
			l7.AutoSize = true;
			l7.Location = new Point (numericUpDown7.Left, numericUpDown7.Bottom);
			l7.Size = numericUpDown7.Size;
			l7.Parent = form1;
			Label l6 = new Label ();
			l6.Text = string.Format("Sampling rate: {0}\nDuration: {1}",sfd,nSec);
			l6.AutoSize = true;
			l6.Location = new Point (10, 470);
			//l6.Size = numericUpDown5.Size;
			l6.Parent = form1;
			NumericUpDown numericUpDown2 = new NumericUpDown ();
			numericUpDown2.Parent = form1;
			numericUpDown2.Maximum = 5m;
			numericUpDown2.Minimum = -5m;
			numericUpDown2.Value = 0m;
			numericUpDown2.Increment = 0.5m;
			numericUpDown2.DecimalPlaces=1;
			numericUpDown2.Location = new Point (10, 120);
			Label l2 = new Label ();
			l2.Text = "Repolarization threshold";
			l2.AutoSize = true;
			l2.Location = new Point (numericUpDown2.Left, numericUpDown2.Bottom);
			l2.Size = numericUpDown2.Size;
			l2.Parent = form1;
			DomainUpDown upDown2 = new DomainUpDown ();
			upDown2.Parent = form1;
			upDown2.Items.Add ("yes");
			upDown2.Items.Add ("no");
			upDown2.SelectedIndex = 1;
			upDown2.Location = new Point (10, 420);
			Label l02 = new Label ();
			l02.Text = "measure autocorrelation?";
			l02.AutoSize = true;
			l02.Location = new Point (upDown2.Left, upDown2.Bottom);
			l02.Size = upDown2.Size;
			l02.Parent = form1;
			if (Indices [0] == 0) {
				DomainUpDown upDown1 = new DomainUpDown ();
				upDown1.Parent = form1;
				upDown1.Items.Add ("yes");
				upDown1.Items.Add ("no");
				upDown1.SelectedIndex = 0;
				upDown1.Location = new Point (10, 170);
				Label l0 = new Label ();
				l0.Text = "Recalibration";
				l0.AutoSize = true;
				l0.Location = new Point (upDown1.Left, upDown1.Bottom);
				l0.Size = upDown1.Size;
				l0.Parent = form1;
				form1.ShowDialog ();
				recalibTrigger = upDown1.SelectedIndex;
			}
			else {
				form1.ShowDialog ();
			}
			// set the parameters
			threshold = (int)(numericUpDown1.Value*AmpScale);
			//thr5 = (threshold + 2) / 2;
			df = (int)numericUpDown3.Value;
			AHPthr = (int)(numericUpDown2.Value*AmpScale);
			CutPre = (int)(numericUpDown4.Value/FrameDt+0.1m);
			CutPost = (int)(numericUpDown5.Value/FrameDt+0.1m);
			tSmth = (int)(numericUpDown7.Value/FrameDt+0.1m);
			tSmth1 = tSmth - 1;
			ACF=(upDown2.SelectedIndex==0);
			form1.Dispose();
			Avgs2= new int[NChannels,NFblocks];
			if (sf / 5000 > 3) {
				NFblocks = sf / 5000;
			}
			FiltNorm = 2*(8 * (NFblocks - 1) + 2);
			HF = (sf > 12000);
			if (HF) {
				TauFiltOffset = TauFilt / 2 + 4 * TauFilt1 *(NFblocks);
			} else {
				TauFiltOffset = TauFilt / 2 + 2 * TauFilt1*(NFblocks+1);
			}
			if (df >= 0) {
				if (df == 0) {
					df += 1;
				}
				dfAbs = df;
				dfSign = 1;
				Slmax = (sf / 1000)/df +1;
				CutOffset = (sf / 1000 + CutPre)/df-tSmth/2;//from where to start cutting...
				//CutAfter = (sf / 1002)/df;
				//CutAfterLong = (sf / 501)/df;
				tCut = (CutPre + 1 + CutPost) / df;//6 +(sf / 1002 + sf / 1000)/df;//cutout length
				tCutLong = tCut + sfi / df;//6 +(sf / 501 + sf / 1000)/df;// long cutout length
				Sln0= (sf / 2334)/df;//3340//835*12/5
				Sampling = sf/df;
				SqIglobal = 0;
				tx = 0;
				ti = df * Math.Max(CutOffset,8);//8 needed for correlation estimate!
				tf =  df * (tInc -Math.Max(tCutLong-CutOffset,Math.Max(8,TauFiltOffset+2*TauFilt1)));
				tm=(tInc-1)*df-tf+ti;
				tdSmth1=-(tSmth-1);
				tdSmth=-tSmth;
				tms=(tInc-tSmth1-1)*df-tf+ti;
				tfiA = (tf - ti) / dfAbs;
				ty=df*(tInc-1);
				t0x = 0;
				dfTI [1] = tf-ti;
			}
			else {
				dfAbs = -df;
				dfSign = -1;
				Slmax = -(sf / 1000)/df+1;
				CutOffset =(sf / 1000 - CutPre)/df+(tSmth-1)/2;// 4 +(sf/835)/df;//should be ok if negative
				//CutAfter = 2 -(sf / 1002 + sf/835 + sf / 1000)/df;
				//CutAfterLong = 2 -(sf / 1002 + sf/835 + sf / 1000)/df;
				tCut=-(CutPre + 1 + CutPost) / df;//6 -(sf / 1002 + sf / 1000)/df;
				tCutLong=-(CutPre + 1 + CutPost) / df;//same as tCut here; long spikes do not make sense
				//tCutLong0=-(CutPre + 1 + CutPost + sfi) / df;//6 -(sf / 501 + sf / 1000)/df;
				Sln0= -(sf / 2334)/df;
				Sampling = -sf/df;
				SqIglobal = 0;
				tx = (-tInc+1)*df;
				ti =  df *(Math.Max(tCutLong-CutOffset,9)-tInc);//backward starts with -1...
				tf = Math.Max(Math.Max(7,TauFiltOffset+2*TauFilt1-1),CutOffset)*dfAbs;//-df * CutOffset;
				tm=ti-tf;
				tdSmth1=(tSmth-1);
				tdSmth=tSmth;
				tms=-(tSmth1)*df-tf+ti;
				tfiA = (tf - ti) / dfAbs;
				ty=0;
				t0x = nFrames + df*tInc;
				dfTI [1] = ti-tf;//positive
			}
			TauFiltOffset *= dfSign;
			Subsample = (dfAbs > 1);
			tQm0 = Sampling / 2000;
			tQmi = (Sampling + 1000) / 2000+tQm0;
			tQmf = (Sampling * 6) / 5000+tQmi;
			//tQml = tQmf + tQmi+tQm0;
			tQm= new int[tQmf];
			tQmLong= new int[tQmf];
			tQmX= new int[tQmf];
			tQmXLong= new int[tQmf];
			for (int i=tQm0; i<tQmi; i++) {
				tQm [i] = -CutOffset + i-tQm0;//-8:-4
				tQmLong [i] = -CutOffset + i-tQm0;
				tQmX [i] = i-tQm0;//-8:-4
				tQmXLong [i] = i-tQm0;
			}
			for (int i=tQmi; i<tQmf; i++) {
				tQm [i] = -CutOffset+tCut-tQmf+ i;//4:12
				tQmLong [i] = -CutOffset+tCutLong-tQmf + i;
				tQmX [i] = tCut +i-tQmf;//-8:-4
				tQmXLong [i] = tCutLong+i-tQmf;
			}
			tQmm = (tQmf) / 2;
			tQmA= new int[tQmf];
			Lspike = Sampling / 5000 + 3;
			Lsw = (Lspike - 1);
			Li = CutPre - CutPre / 3 - Lspike / 2;//-4:4 how likely Overlap?
			Lf = CutPre + CutPre / 3 - Lspike / 2 +1;
			//Ll = Lf - Li + Lspike - 1;
			//LiX = - CutPre / 3 - Lspike / 2;//-4:4 how likely Overlap?
			//LfX = CutPre / 3 - Lspike / 2 +1;
			Lw=new int[Lspike];
			Lw [0] = 1;
			Lw [Lspike - 1] = 1;
			for (int i =1; i<Lspike-1; i++) {
				Lw [i] = 2;
			}
			Lmx= 0;
			tShape= new int[tCutLong];
			dfTI [0] = df;
			dfTI [2] = tInc*dfAbs;//positive
			Ascale = (Ascale0 / tSmth / tSmth1/2) *2* tSmth * tSmth1;
			AscaleV = Ascale / tSmth;
			AscaleG = Ascale / tSmth1/2;
			AglobalSmth=new int[tInc];
			for (int i=0; i<tInc; i++) {
				AglobalSmth [i] = 2047*Ascale;
			}
			vSmth=new int[NChannels];
			//TauF1 = 2 * TauFilt1 * dfAbs;
			// open a box that logs acticity during detection
			logwindow = new Form();
			logwindow.Size = new Size(210, 180);
			logbox = new TextBox();
        	logbox.Parent = logwindow;
        	logbox.Dock = DockStyle.Fill;
        	logbox.Multiline = true;
			logwindow.Visible = true;
			logbox.Visible = true;
			logbox.AppendText("Detecting spikes...\n");
			logbox.Update();
			logwindow.Update();


			for (int i=0;i<NChannels;i++){
		   		Qd[i]=600;
				Qm[i]=0;
				for (int ij=0; ij<dtMx; ij++) {
					Qdiff [i, ij] = 0;
					Qmax[i,ij]=0;
				}
				QmaxE[i]=0;
				SqIv[i]=0;//sum of squared channel increments
				//SIp[i]=0;
				for (int ii=0; ii<13; ii++) {
					SIprod[i,ii]=0;//sum of product of global and channel voltage increments
				}
				for (int iii=0; iii<3; iii++) {
					Avgs1 [i, iii] = 4094;//sum of product of global and channel voltage increments
				}
				for (int iii=0; iii<NFblocks; iii++) {
					Avgs2 [i, iii] = 4094;//sum of product of global and channel voltage increments
				}
				QmPreD [i] = 0;
				Avgs3 [i] = 2047*Ascale;
				FVbias[i]=0;
				Vbias[i]=0;
				FVsbias[i]=0;
				Vsqbias[i]=0;
				A[i]=0;
				Sl4[i]=0;
				Sl4x[i]=false;
				Z4next[i]=0;
				Z5next[i]=0;
				AHP4[i]=false;
				Amp4[i]=0;
				Sl5[i]=0;
				Sl5x[i]=false;
				AHP5[i]=false;
				Amp5[i]=0;
			}
			return dfTI;
		}


		
		public void openSpikeFile (string name) {
			fs = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			w = new StreamWriter(fs);
		}
		
		public void openShapeFile (string name) {
			fsShapes = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wShapes = new StreamWriter(fsShapes);
		}
		
		public void openSpikeXFile (string name) {
			fsX = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wX = new StreamWriter(fsX);
		}
		
		public void openShapeXFile (string name) {
			fsShapesX = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wShapesX = new StreamWriter(fsShapesX);
		}

		public void openInfoFile (string name) {
			fsInfo = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wInfo = new StreamWriter(fsInfo);
		}

		public void openMeanFile (string name) {
			fsMean = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wMean = new StreamWriter(fsMean);
		}

		public void AvgVoltageDefault (short[][] vm, long t0, int t)//want to compute an approximate 33 percentile
		{
			//can average over 2 consecutive frames
			//each time called, I should take the next 4 frames of vm (all channels)
			//would need to correct for Aglobal
			//iterates over 8 frames
			int P1;
			int P2;
			int Plow = 0;//lowest value to discard
			int Plowx = 0;//temporary variable
			int k;
			int kk;
			bool Px;
			if (Subsample) {
				if (HF) {
					k = (int)((t + t0) / (8 * TauFilt1)) % NFblocks;
					kk = (int)((t0 + t) / (TauFilt1 * 2)) % 4;
				} else {
					k = (int)((t + t0) / (TauFilt1 * 8)) % NFblocks;
					if ((((t0 + t) / (TauFilt1 * 2)) % 2) != 0) {
						kk = 2;
					} else {
						kk = (int)(((t0 + t) / (4 * TauFilt1)) % 2);
					}
				}
				int knext = (k + (NFblocks - 1)) % NFblocks;
				int klast = (k + 1) % NFblocks;
				for (int i=1; i<NChannels; i++) {//loop across channels
					//want something that works in most cases
					Px = false;
					//assume that I already have the right value
					P1 = 2*vm [i] [t]-Aglobal[t/dfAbs];
					for (int tt=1; tt<TauFilt1; tt++) {//this function wastes most of the time
						Plowx = 2*vm [i] [t + tt * df]-Aglobal[(t+tt*df)/dfAbs];//factor of 3!
						if (Plowx < P1) {
							if (!Px) {
								Px = true;
								Plow = Plowx;
							} else {
								P1 = Plow;
								Plow = Plowx;
							}
						}
					}
					if (!Px) {//P1 was the lowest value
						P1 = Plowx;
						for (int tt=TauFilt1-2; tt>0; tt--) {
							Plowx = 2*vm [i] [t + tt * df]-Aglobal[(t+tt*df)/dfAbs];
							if (Plowx < P1) {
								P1 = Plowx;
							}
						}
					}
					//outliers
					if (((P1 + 13000) % 16000) < 10000) {
						P1 = 0;
					}
					//how about average over four frames, take weighted average (first and last weighted by 1/2 or so)
					//i.e. for intermediate values, just change the weighting
					//(advantage: less memory wastage and only need to compute median every fourth frame)
					//have one block of four values to estimate (make area to average shorter), take weighted average of first and last block
					if (kk < 2) {
						Avgs3 [i] += (Avgs2 [i, knext] - Avgs2 [i, klast]) * Ascale / FiltNorm;
						Avgs1 [i, kk] = 2*P1;
					} else if (HF & (kk<3)) {
						Avgs3 [i] += (Avgs2 [i, knext] - Avgs2 [i, klast]) * Ascale / FiltNorm;
						Avgs1 [i, kk] = 2*P1;
					} else {
						//assume that I already have the right value
						Px = false;
						P2 = 2*P1;
						for (int tt=0; tt<3; tt++) {//this function wastes most of the time
							Plowx = Avgs1 [i, tt];
							if (Plowx < P2) {
								if (!Px) {
									Px = true;
									Plow = Plowx;
								} else {
									P2 = Plow;
									Plow = Plowx;
								}
							}
						}
						if (!Px) {//P1 was the lowest value
							P2 = Plowx;
							for (int tt=2; tt>=0; tt--) {
								Plowx = Avgs1 [i, tt];
								if (Plowx < P2) {
									P2 = Plowx;
								}
							}
						}
						Avgs2 [i, klast] = P2;
						if (!HF) {
							Avgs1 [i, 2] = 2*P1;
						}
						//to avoid accumulating numerical errors (not sure whether necessary)
						if (NFblocks == 3) {
							Avgs3 [i] = (P2 + 4 * (Avgs2 [i, knext] + Avgs2 [i, k])) * Ascale/ FiltNorm;
						} else {
							Avgs3 [i] = -3*P2;
							for (int tt=0; tt<NFblocks; tt++) {
								Avgs3 [i] += 4 * Avgs2 [i, tt];
							}
							Avgs3 [i] *= Ascale;
							Avgs3 [i] /=FiltNorm;
						}
					}
				}
			} else {
				if (df > 0) {
					if (HF) {
						k = (int)((t + t0) / (8 * TauFilt1)) % NFblocks;
						kk = (int)((t0 + t) / (2 * TauFilt1)) % 4;
					} else {
						k = (int)((t + t0) / (4 * TauFilt1)) % NFblocks;//overlapping blocks
						if ((((t0 + t) / (2 * TauFilt1)) % 2) != 0) {
							kk = 2;
						} else {
							kk = (int)(((t0 + t) / (4 * TauFilt1)) % 2);
						}
					}
					int knext = (k + (NFblocks - 1)) % NFblocks;
					int klast = (k + 1) % NFblocks;
					//Console.WriteLine ("{0} {1} {2}", k, knext, klast);
					//int FiltNorm = 18;
					//int Pl = 0;//lowest value to discard
					//int Plx = 0;//temporary variable
					//int[] Px1= new int[TauFilt1];
					for (int i=1; i<NChannels; i++) {//loop across channels
						//want something that works in most cases
						Px = false;
						//assume that I already have the right value
						P1 = 2*(vm [i] [t] + vm [i] [t + 1])-Aglobal[t]-Aglobal[t+1];
						for (int tt=1; tt<TauFilt1; tt++) {//this function wastes most of the time
							Plowx = 2*(vm [i] [t + 2 * tt] + vm [i] [t + 2 * tt + 1])-Aglobal[t+2*tt]-Aglobal[t+2*tt+1];
							if (Plowx < P1) {
								if (!Px) {
									Px = true;
									Plow = Plowx;
								} else {
									P1 = Plow;
									Plow = Plowx;
								}
							}
						}
						if (!Px) {//P1 was the lowest value
							P1 = Plowx;
							for (int tt=TauFilt1-2; tt>0; tt--) {
								Plowx = 2*(vm [i] [t + 2 * tt] + vm [i] [t + 2 * tt + 1])-Aglobal[t+2*tt]-Aglobal[t+2*tt+1];
								if (Plowx < P1) {
									P1 = Plowx;
								}
							}
						}
						//outliers
						if (((P1 + 26000) % 32000) < 20000) {
							P1 = 0;
						}
						//how about average over four frames, take weighted average (first and last weighted by 1/2 or so)
						//i.e. for intermediate values, just change the weighting
						//(advantage: less memory wastage and only need to compute median every fourth frame)
						//have one block of four values to estimate (make area to average shorter), take weighted average of first and last block
						//Avgs1 [i, kk] = P1;
						if (kk < 2) {
							Avgs3 [i] += (Avgs2 [i, knext] - Avgs2 [i, klast]) * Ascale / FiltNorm;
							Avgs1 [i, kk] = P1;
						} else if (HF & (kk < 3)) {
							Avgs3 [i] += (Avgs2 [i, knext] - Avgs2 [i, klast]) * Ascale / FiltNorm;
							Avgs1 [i, kk] = P1;
						} else {
							//assume that I already have the right value
							P2 = P1;
							Px = false;
							for (int tt=0; tt<3; tt++) {//this function wastes most of the time
								Plowx = Avgs1 [i, tt];
								if (Plowx < P2) {
									if (!Px) {
										Px = true;
										Plow = Plowx;
									} else {
										P2 = Plow;
										Plow = Plowx;
									}
								}
							}
							if (!Px) {//P1 was the lowest value
								P2 = Plowx;
								for (int tt=1; tt>=0; tt--) {
									Plowx = Avgs1 [i, tt];
									if (Plowx < P2) {
										P2 = Plowx;
									}
								}
							}
							Avgs2 [i, klast] = P2;//klast will be knext after one iteration
							if (!HF) {
								Avgs1 [i, 2] = P1;
							}
							//to avoid accumulating numerical errors (not sure whether necessary)
							if (NFblocks == 3) {
								Avgs3 [i] = (P2 + 4 * (Avgs2 [i, knext] + Avgs2 [i, k])) * Ascale / FiltNorm;
							} else {
								Avgs3 [i] = -3 * P2;
								for (int tt=0; tt<NFblocks; tt++) {
									Avgs3 [i] += 4 * Avgs2 [i, tt];
								}
								Avgs3 [i] *= Ascale;
								Avgs3 [i] /= FiltNorm;
							}
						}
					}
				} else {
					if (HF) {
						k = (int)((t + t0) / (8 * TauFilt1)) % NFblocks;
						kk = (int)((t0 + t) / (2 * TauFilt1)) % 4;
					} else {
						k = (int)((t + t0) / (4 * TauFilt1)) % NFblocks;//overlapping blocks
						if ((((t0 + t) / (2 * TauFilt1)) % 2) != 0) {
							kk = 2;
						} else {
							kk = (int)(((t0 + t) / (4 * TauFilt1)) % 2);
						}
					}
					int knext = (k + (NFblocks - 1)) % NFblocks;
					int klast = (k + 1) % NFblocks;
					//Console.WriteLine ("{0} {1} {2}", k, knext, klast);
					//int FiltNorm = 18;
					//int Pl = 0;//lowest value to discard
					//int Plx = 0;//temporary variable
					//int[] Px1= new int[TauFilt1];
					for (int i=1; i<NChannels; i++) {//loop across channels
						//want something that works in most cases
						Px = false;
						//assume that I already have the right value
						P1 = 2*(vm [i] [t] + vm [i] [t - 1])-Aglobal[t]-Aglobal[t-1];
						for (int tt=1; tt<TauFilt1; tt++) {//this function wastes most of the time
							Plowx = 2*(vm [i] [t - 2 * tt] + vm [i] [t - 2 * tt - 1])-Aglobal[t-2*tt]-Aglobal[t-2*tt-1];
							if (Plowx < P1) {
								if (!Px) {
									Px = true;
									Plow = Plowx;
								} else {
									P1 = Plow;
									Plow = Plowx;
								}
							}
						}
						if (!Px) {//P1 was the lowest value
							P1 = Plowx;
							for (int tt=TauFilt1-2; tt>0; tt--) {
								Plowx = 2*(vm [i] [t - 2 * tt] + vm [i] [t - 2 * tt - 1])-Aglobal[t-2*tt]-Aglobal[t-2*tt-1];
								if (Plowx < P1) {
									P1 = Plowx;
								}
							}
						}
						//outliers
						if (((P1 + 26000) % 32000) < 20000) {
							P1 = 0;
						}
						//how about average over four frames, take weighted average (first and last weighted by 1/2 or so)
						//i.e. for intermediate values, just change the weighting
						//(advantage: less memory wastage and only need to compute median every fourth frame)
						//have one block of four values to estimate (make area to average shorter), take weighted average of first and last block
						//Avgs1 [i, kk] = P1;
						if (kk < 2) {
							Avgs3 [i] += (Avgs2 [i, knext] - Avgs2 [i, klast]) * Ascale / FiltNorm;
							Avgs1 [i, kk] = P1;
						} else if (HF & (kk < 3)) {
							Avgs3 [i] += (Avgs2 [i, knext] - Avgs2 [i, klast]) * Ascale / FiltNorm;
							Avgs1 [i, kk] = P1;
						} else {
							//assume that I already have the right value
							P2 = P1;
							Px = false;
							for (int tt=0; tt<3; tt++) {//this function wastes most of the time
								Plowx = Avgs1 [i, tt];
								if (Plowx < P2) {
									if (!Px) {
										Px = true;
										Plow = Plowx;
									} else {
										P2 = Plow;
										Plow = Plowx;
									}
								}
							}
							if (!Px) {//P1 was the lowest value
								P2 = Plowx;
								for (int tt=1; tt>=0; tt--) {
									Plowx = Avgs1 [i, tt];
									if (Plowx < P2) {
										P2 = Plowx;
									}
								}
							}
							Avgs2 [i, klast] = P2;//klast will be knext after one iteration
							if (!HF) {
								Avgs1 [i, 2] = P1;
							}
							//to avoid accumulating numerical errors (not sure whether necessary)
							if (NFblocks == 3) {
								Avgs3 [i] = (P2 + 4 * (Avgs2 [i, knext] + Avgs2 [i, k])) * Ascale / FiltNorm;
							} else {
								Avgs3 [i] = -3 * P2;
								for (int tt=0; tt<NFblocks; tt++) {
									Avgs3 [i] += 4 * Avgs2 [i, tt];
								}
								Avgs3 [i] *= Ascale;
								Avgs3 [i] /= FiltNorm;
							}
						}
					}
				}
			}
			//Console.WriteLine ("{0} {1} {2}",Avgs2[2120,0],Avgs2[2120,1], Avgs2[2120,2]);
		}

		//should use the slow filter here as well, to avoid transients
		public void InitialEstimation (short[][] vm, long t0) {//use this to get a better initial estimate of Qd. only fast transients.
			int tA;
			if (t0 == t0x) {
				//estimate Aglobal
				for (int t=tx; dfSign*t<dfSign*ty; t+=df) {//loop over data, will be removed for an online algorithm
					tA = t / dfAbs;
					for (int i=1; i<NChannels; i++) {//loop across channels
						Slice [i] = (vm [i] [t]) % 4095 + (vm [i] [t + df]) % 4095;
					}
					Array.Sort (Slice);
					Aglobal [tA] = Slice [NChannels / 2];
				}
				for (int t=tx/dfAbs+dfSign; dfSign*t<ty/df; t+=dfSign) {
					Aglobaldiff [t] = Aglobal [t] - Aglobal [t - dfSign];
				}
				for (int t=tx/dfAbs; dfSign*t<ty/df-dfSign*tSmth1; t+=dfSign) {
					//tA = t / dfAbs;
					AglobalSmth [t] = Aglobal [t];
					for (int ii=1; ii<tSmth1; ii++) {
						AglobalSmth [t] += Aglobal [t + ii*dfSign];
					}
					AglobalSmth [t] *= AscaleG;
				}
				for (int t=tx/dfAbs+dfSign; dfSign*t<ty/df-dfSign*tSmth1; t+=dfSign) {
					AglobalSdiff [t] = AglobalSmth [t] - AglobalSmth [t - dfSign];
				}
				//initialize slow filter
				if (HF) {
					for (int t=tx; dfSign*t<dfSign*tx+NFblocks*4*TauFilt1*2; t+=2*TauFilt1*dfSign) {
						AvgVoltageDefault (vm, t0, t);//-4+4*dfSign
					}
				} else {
					for (int t=tx; dfSign*t<dfSign*tx+(NFblocks+1)*4*TauFilt1; t+=2*TauFilt1*dfSign) {
						AvgVoltageDefault (vm, t0, t);
					}
				}
				for (int t=tx; dfSign*t<dfSign*ti; t+=df) {
					tA = t / dfAbs;
					for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
						vSmth [i] = vm [i] [t];
						for (int ii=1; ii<tSmth; ii++) {
							vSmth [i] += vm [i] [t + ii * df];
						}
						vSmth [i] *= AscaleV;
						//CHANNEL OUT OF LINEAR REGIME
						if (((vm [i] [t + 2*df] + 4) % 4096) < 10) {
							if (A [i] < artT) {//reset only when it starts leaving the linear regime
								A [i] = artT;
							}
						} else {
							Qm [i] = (2 * (Qm [i] + Qd [i]) + vSmth [i] - AglobalSmth[tA]) / 3;//update Qm
						}
					}
				}
				//shift Aglobal entries
				for (int t=tx/dfAbs; dfSign*t<tm/df; t+=dfSign) {
					Aglobal [t + tfiA] = Aglobal [t];
					Aglobaldiff [t + tfiA] = Aglobaldiff [t];
				}
				for (int t=tx/dfAbs; dfSign*t<tms/df; t+=dfSign) {
					AglobalSmth [t + tfiA] = AglobalSmth [t];
					AglobalSdiff [t + tfiA] = AglobalSdiff [t];
				}
			}
			//shift Aglobal entries
			for (int t=tx/dfAbs; dfSign*t<tm/df; t+=dfSign) {
				Aglobal [t] = Aglobal [t + tfiA];
				Aglobaldiff [t] = Aglobaldiff [t + tfiA];
			}
			for (int t=tx/dfAbs; dfSign*t<tms/df; t+=dfSign) {
				AglobalSmth [t] = AglobalSmth [t + tfiA];
				AglobalSdiff [t] = AglobalSdiff [t + tfiA];
			}
			//new Aglobal entries
			for (int t=tm; dfSign*t<dfSign*ty; t+=df) {
				tA = t / dfAbs;
				for (int i=1; i<NChannels; i++) {//loop across channels
					Slice[i]=(vm [i][t])%4095+(vm[i][t+df])%4095;
				}
				Array.Sort(Slice);
				Aglobal[tA]=Slice[NChannels / 2];
				Aglobaldiff [tA] = Aglobal [tA] - Aglobal [tA - dfSign];
				AglobalSmth [tA+tdSmth1] = Aglobal [tA+tdSmth1];
				for (int ii=1; ii<tSmth1; ii++) {
					AglobalSmth [tA+tdSmth1] += Aglobal [tA+tdSmth1 + ii*dfSign];
				}
				AglobalSmth [tA+tdSmth1] *= AscaleG;
				AglobalSdiff[tA+tdSmth1]=AglobalSmth[tA+tdSmth1]-AglobalSmth[t+tdSmth];
			}
			//avoid cumulation of numerical errors
			for (int i=1; i<NChannels; i++) {//loop across channels
				vSmth [i] = vm [i] [ti - df];
				for (int ii=0; ii<tSmth1; ii++) {
					vSmth [i] += vm [i] [ti + ii * df];
				}
				vSmth [i] *= AscaleV;
			}
			for (int t=ti; dfSign*t<dfSign*tf; t+=df) {//loop over data, will be removed for an online algorithm
				dt = (dt+1) % dtMx;
				tA = t / dfAbs;
				if ((t0+t)%(2*TauFilt1)==0) {
					//update Avgs3
					AvgVoltageDefault (vm, t0, t+TauFiltOffset);
				}
				for (int i=1; i<NChannels; i++) {//update vSmth
					vSmth [i] += (vm [i] [t+ tSmth1 * df]-vm [i] [t - df]) *AscaleV;
					QmPre[i]=vSmth[i]-AglobalSmth[tA];
				}
				// SPIKE DETECTION
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					//CHANNEL OUT OF LINEAR REGIME
					if (((vm [i] [t + df] + 4) % 4096) < 10) {
						if (A [i] < artT) {//reset only when it starts leaving the linear regime
							A [i] = artT;
						}
					}
					//DEFAULT OPERATIONS
					else if (A [i] == 0) {
						QmPreD [i] += (Avgs3[i] - QmPre [i]-QmPreD [i]) /(TauFilt);
						Vbias[i]=FVbias[i]*AglobalSdiff[tA]/Sampling;
						//vm[i][t-df]+vm[i][t]+vm[i][t+df]
						Qdiff [i, dt]= (QmPre [i]+QmPreD [i])-Qm[i]-Vbias[i];//difference between ADC counts and Qm
						if ((AglobalSdiff[tA] * (Qdiff [i, dt] - Qdiff [i, dtPre])) != 0) {
							if ((AglobalSdiff[tA] > 0) == ((Qdiff [i, dt] - Qdiff [i, dtPre]) > 0)) {//(SIp[i]>0) {
								FVbias [i]++;
							} else {
								FVbias [i]--;
							}//Qdiff negative-->Ascale!!!;
						}
						//Qdiff [i, dt] = (vm [i] [t - df] + vm [i] [t] + vm [i] [t + df] - Aglobal[tA]) * Ascale - Qm [i];//difference between ADC counts and Qm
						//UPDATE Qm and Qd
						if (Qdiff [i, dt] > 0) {
							if (Qdiff [i, dt] > Qd [i]) {
								Qm [i] += Qd [i] / Tau_m0;
								if (Qdiff [i, dt] < (5 * Qd [i])) {
									Qd [i]++;
								} else if ((Qd [i] > Qdmin) & (Qdiff [i, dt] > (6 * Qd [i]))) {
									Qd [i]--;
								}
							} else if (Qd [i] > Qdmin) {//set a minimum level for Qd
								Qd [i]--;
							}
						} else if (Qdiff [i, dt] < -Qd [i]) {
							Qm [i] -= Qd [i] / Tau_m0 / 2;

						}
					}
					//AFTER CHANNEL WAS OUT OF LINEAR REGIME
					else {
						Qm [i] = (2 * (Qm [i] + Qd [i]) + vSmth [i] - AglobalSmth[tA]) / 3;//update Qm
						A [i]--;
					}
				}
			}
			if (df > 0) {
				if (t0 >=199*(dfTI[2] - dfTI[1])) {
					for (int i=0; i<NChannels; i++) {
						Qm [i] = 0;
						for (int ij=0; ij<dtMx; ij++) {
							Qdiff [i, ij] = 0;
						}
						A [i] = 0;
					}
				}
			} else {
				if (t0 <=t0x-199*(dfTI[2] - dfTI[1])) {
					for (int i=0; i<NChannels; i++) {
						Qm [i] = 0;
						for (int ij=0; ij<dtMx; ij++) {
							Qdiff [i, ij] = 0;
						}
						A [i] = 0;
					}
				}
			}
		}

		public void StartDetection (short[][] vm, long t0, long nFrames, double nSec, double sfd, int[] Indices){
			w.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wShapes.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wX.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wShapesX.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wInfo.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wMean.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			//write some info
			wInfo.WriteLine("# Number of frames:\n{0}", nFrames/dfAbs);
			wInfo.WriteLine("# Duration (s):\n{0}", nSec);
			wInfo.WriteLine("# Sampling rate:\n{0}", sfd/dfAbs);
			wInfo.WriteLine("# Threshold scaling:\n{0}", AmpScale);
			wInfo.WriteLine("# Amplitude scaling:\n{0}", Ascale);
			wInfo.WriteLine("# Detection threshold*{0}:\n{1}", AmpScale, threshold);
			wInfo.WriteLine("# Repolarization threshold*{0}:\n{1}", AmpScale, AHPthr);
			wInfo.WriteLine("# Recalibration trigger:\n{0}", recalibTrigger);
			wInfo.WriteLine("# Cutouts:\n{0} {1} {2} {3} {4}", CutPre, CutPost, tCut, tCutLong, df);
			wInfo.WriteLine("# Smoothing window (detection):\n{0}", tSmth);
			wInfo.WriteLine("# Smoothing window (amplitudes):\n{0}", Lspike);
			wInfo.WriteLine ("# Recording channels:");
			for (int i=0; i<Indices.Length; i++) {
				wInfo.WriteLine ("{0}", Indices [i]);
			}
			wInfo.WriteLine ("# Recording channels4:");
			for (int i=0; i<NChannels; i++) {
				for (int j=0; j<4;j++){
					wInfo.Write("{0} ", ChInd4a[i][j]);
				}
				for (int j=0; j<8;j++){
					wInfo.Write("{0} ", ChInd4b[i][j]);
				}
				wInfo.WriteLine ();
			}
			wInfo.WriteLine ("# Recording channels5:");
			for (int i=0; i<NChannels; i++) {
				wInfo.Write("{0} ", ChInd5[i]);
				for (int j=0; j<4;j++){
					wInfo.Write("{0} ", ChInd5a[i][j]);
				}
				for (int j=0; j<4;j++){
					wInfo.Write("{0} ", ChInd5b[i][j]);
				}
				wInfo.WriteLine ();
			}
			wInfo.WriteLine("# Recalibration events:");
			int tA;
			//estimate Aglobal
			for (int t=tx; dfSign*t<dfSign*ty; t+=df) {//loop over data, will be removed for an online algorithm
				tA = t / dfAbs;
				for (int i=1; i<NChannels; i++) {//loop across channels
					Slice [i] = (vm [i] [t]) % 4095 + (vm [i] [t+df]) % 4095;
				}
				Array.Sort (Slice);
				Aglobal[tA] = Slice [NChannels / 2];
			}
			for (int t=tx/dfAbs+dfSign; dfSign*t<ty/df; t+=dfSign) {
				Aglobaldiff [t] = Aglobal [t] - Aglobal [t - dfSign];
			}
			for (int t=tx/dfAbs; dfSign*t<ty/df-dfSign*tSmth1; t+=dfSign) {
				//tA = t / dfAbs;
				AglobalSmth [t] = Aglobal [t];
				for (int ii=1; ii<tSmth1; ii++) {
					AglobalSmth [t] += Aglobal [t + ii*dfSign];
				}
				AglobalSmth [t] *= AscaleG;
			}
			for (int t=tx/dfAbs+dfSign; dfSign*t<ty/df-dfSign*tSmth1; t+=dfSign) {
				AglobalSdiff [t] = AglobalSmth [t] - AglobalSmth [t - dfSign];
			}
			//initialize slow filter
			if (HF) {
				for (int t=tx; dfSign*t<dfSign*tx+NFblocks*4*TauFilt1*2; t+=2*TauFilt1*dfSign) {
					AvgVoltageDefault (vm, t0, t);//-4+4*dfSign
				}
			} else {
				for (int t=tx; dfSign*t<dfSign*tx+(NFblocks+1)*4*TauFilt1; t+=2*TauFilt1*dfSign) {
					AvgVoltageDefault (vm, t0, t);
				}
			}
			for (int t=tx; dfSign*t<dfSign*ti; t+=df) {
				tA = t / dfAbs;
				//Console.Write("{0} ", t+t0);
				dtPre = dt;
				dt = (dt+1) % dtMx;
				dtEx = dtE;
				dtE = (dtE + 1) % dtEMx;
				if ((t0 + t) % (2 * TauFilt1) == 0) {
					AvgVoltageDefault (vm, t0, t + TauFiltOffset);
				}
				//now have to change voltage
				//can use Qm?
				for (int i=1; i<NChannels; i++) {//loop across channels
					vSmth [i] = vm [i] [t];
					for (int ii=1; ii<tSmth; ii++) {
						vSmth [i] += vm [i] [t + ii * df];
					}
					vSmth [i] *= AscaleV;
					QmPre[i]=vSmth[i]-AglobalSmth[tA];
					//QmPre[i]=(vm[i][t]+vm [i][t+df]+vm[i][t+2*df]-Aglobal[tA])*Ascale;
				//	//QmPreD [i] += (Avgs3[i] - QmPre [i]-QmPreD [i]) /(TauFilt);
				//	Slice [i] = QmPre [i]+QmPreD [i];//((vm[i][t])%4095+(vm[i][t+df])%4095+(vm[i][t+2*df])%4095);
				}

				//Array.Sort(Slice);
				//Aglobaldiff=Slice[NChannels / 2]-Aglobal;
				//Aglobal=Slice[NChannels / 2];
				wMean.WriteLine("{0}", Aglobal[tA]);
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					//CHANNEL OUT OF LINEAR REGIME
					if (((vm[i][t+2*df]+4)%4096)<10) {
						if (A[i]<artT) {//reset only when it starts leaving the linear regime
							A[i]=artT;
							QmPreD[i] = 0;
						}
					}
					else {
						QmPreD[i]  += (Avgs3[i] - QmPre [i]-QmPreD [i]) /(TauFilt);
						Qm[i]=(2*(Qm[i]+Qd[i])+QmPre [i]+QmPreD [i])/3;//update Qm vm[i][t]+vm[i][t+df]+vm[i][t+2*df]
					}
				}
			}
			//shift Aglobal entries
			for (int t=tx/dfAbs; dfSign*t<tm/df; t+=dfSign) {
				Aglobal [t + tfiA] = Aglobal [t];
				Aglobaldiff [t + tfiA] = Aglobaldiff [t];
			}
			for (int t=tx/dfAbs; dfSign*t<tms/df; t+=dfSign) {
				AglobalSmth [t + tfiA] = AglobalSmth [t];
				AglobalSdiff [t + tfiA] = AglobalSdiff [t];
			}
		}

		public void skipLastReverse (int skipLast) {
			if (df < 0) {
				Console.WriteLine ("{0}", skipLast);
				//ti -= skipLast;
				tx -= skipLast;
				tm -= skipLast;
				tms -= skipLast;
				ti -= skipLast;
				skipLast = 0;
			} else {
				tf -= skipLast;
				ty -= skipLast;
			}
		}

		public void Iterate (short[][] vm, long t0){
			//int qq;
			int a4;//to buffer the difference between ADC counts and Qm
			int a5;//to buffer the difference between ADC counts and Qm
			int b;
			//int[] ChS4 = new int[4];
			//int[] ChS5 = new int[4];
			int ChS4Min;
			int ChS5Min;
			int tA;
			//shift Aglobal entries
			for (int t=tx/ dfAbs; dfSign*t<tm/df; t+=dfSign) {
				Aglobal [t] = Aglobal [t + tfiA];
				Aglobaldiff [t] = Aglobaldiff [t + tfiA];
			}
			for (int t=tx/ dfAbs; dfSign*t<tms/df; t+=dfSign) {
				AglobalSmth [t] = AglobalSmth [t + tfiA];
				AglobalSdiff [t] = AglobalSdiff [t + tfiA];
			}
			//new Aglobal entries
			for (int t=tm; dfSign*t<dfSign*ty; t+=df) {
				tA = t / dfAbs;
				for (int i=1; i<NChannels; i++) {//loop across channels
					Slice[i]=(vm [i][t])%4095+(vm[i][t+df])%4095;
				}
				Array.Sort(Slice);
				Aglobal[tA]=Slice[NChannels / 2];
				Aglobaldiff [tA] = Aglobal [tA] - Aglobal [tA - dfSign];
				AglobalSmth [tA+tdSmth1] = Aglobal [tA+tdSmth1];
				for (int ii=1; ii<tSmth1; ii++) {
					AglobalSmth [tA+tdSmth1] += Aglobal [tA+tdSmth1 + ii*dfSign];
				}
				AglobalSmth [tA+tdSmth1] *= AscaleG;
				AglobalSdiff[tA+tdSmth1]=AglobalSmth[tA+tdSmth1]-AglobalSmth[t+tdSmth];
			}
			//avoid cumulation of numerical errors
			for (int i=1; i<NChannels; i++) {//loop across channels
				vSmth [i] = vm [i] [ti - df];
				for (int ii=0; ii<tSmth1; ii++) {
					vSmth [i] += vm [i] [ti + ii * df];
				}
				vSmth [i] *= AscaleV;
			}
			for (int t=ti; dfSign*t<dfSign*tf; t+=df) {//loop over data, will be removed for an online algorithm
				//Console.Write("{0} ", t+t0);
				dtPre = dt;
				dt=(dt+1)%dtMx;
				dtEx=dtE;
				dtE=(dtE+1)%dtEMx;
				tA = t / dfAbs;
				if ((t0+t)%(2*TauFilt1)==0) {
					//update Avgs3
					AvgVoltageDefault (vm, t0, t+TauFiltOffset);
				}
				for (int i=1; i<NChannels; i++) {//update vSmth
					vSmth [i] += (vm [i] [t+ tSmth1 * df]-vm [i] [t - df]) *AscaleV;
					QmPre[i]=vSmth[i]-AglobalSmth[tA];
				}
				////for (int i=1; i<NChannels; i++) {//loop across channels
				////	QmPre[i]=(vm [i][t]+vm[i][t-df]+vm[i][t+df]-Aglobal[tA])*Ascale;
				////}
				//for (int i=1; i<NChannels; i++) {//loop across channels
				//	QmPre[i]=((vm[i][t-df])%4095+(vm[i][t])%4095+(vm[i][t+df])%4095)*Ascale;
				//	//QmPreD [i] += (Avgs3[i] - QmPre [i]-QmPreD [i]) /(TauFilt);
				//	//if (i == 2120) {
				//	//	Console.WriteLine (QmPreD [i]);
				//	//}
				//	Slice [i] = QmPre [i]+QmPreD [i];//((vm[i][t])%4095+(vm[i][t+df])%4095+(vm[i][t+2*df])%4095);
				//	//Slice[i]=((vm[i][t-df])%4095+(vm[i][t])%4095+(vm[i][t+df])%4095);
				//}
				//Array.Sort(Slice);
				//if (ACF) {
				//	Aglobaldiffold = Aglobaldiff;
				//}
				//Aglobaldiff=Aglobal[tA+dfSign]-Aglobal[tA-dfSign];
				//Aglobal=Slice[NChannels / 2];
				//SqIglobal+=Aglobaldiff[tA]*Aglobaldiff[tA];
				wMean.WriteLine("{0}", Aglobal[tA]);
				// RECALIBRATION EVENTS
				if (recalibTrigger==0){
					if (vm[0][t]<2500) {//write spikes after each recalibration event_newfiles:<2500_oldfiles:<1500
						if (Acal > 2000) {
							wInfo.Write("{0}", (t+t0)/dfAbs);//write time of recalibration event
							for (int i=0; i<NChannels; i++) {//loop across channels
								if (A[i]==0) {
									wInfo.Write(" {0}", Qd[i]);//write variance of recalibration event
								}
								else {
									wInfo.Write (" 0");
								}
							}
							wInfo.WriteLine();
							Console.WriteLine ("{0} sec", (t+t0)/Sampling/dfAbs);// to monitor progress of spike detection
							logbox.AppendText(String.Format("{0} sec\n", (t+t0)/dfAbs/Sampling));
							logwindow.Update();
							logbox.Update();
							Acal=0;//to remember last recalibration event
							for (int i=0; i<NChannels; i++) {
								FVsbias [i] += FVbias [i];//want to see whether this matches the correlation structure
								Vsqbias [i] += (FVbias [i] / 100) * (FVbias [i] / 100);
								SqIglobal+=AglobalSdiff[tA]*AglobalSdiff[tA];
							}
						}
					}
					Acal++;
				}
				else if ((t0+t)%Sampling==0) {//write spikes after every second
					wInfo.Write("{0}", (t+t0)/dfAbs);//write time of recalibration event
					for (int i=0; i<NChannels; i++) {//loop across channels
						if (A[i]==0) {
							wInfo.Write(" {0}", Qd[i]);//write variance of recalibration event
						}
						else {
							wInfo.Write (" 0");
						}
					}
					wInfo.WriteLine();
					Console.WriteLine ("{0} sec", (t+t0)/Sampling/dfAbs);// to monitor progress of spike detection
					logbox.AppendText(String.Format("{0} sec\n", (t+t0)/Sampling/dfAbs));
					logwindow.Update();
					logbox.Update();
				}
				/*
				//for testing
				for (int i=0; i<4096; i++) {
					Qdiff[i,dt]=0;
				}
				for (int i=0; i<24;i++){
					for (int kkk=0; kkk<9; kkk++){
						//Console.WriteLine ("{0} {1} {2} {3}", (t+t0)/Sampling,Iseq[(t+t0-i+256)%256]+Inn[kkk],((t+t0-i+256)%256)/16,(t+t0-i+30)%15);
						Qdiff[Iseq[(t+t0-i+256)%256]+Inn[kkk],dt]=(Iweights[((t+t0-i+256)%256)/16,kkk]*Template[(t+t0-i+256)%16,i])/Ann[(t+t0-i+30)%15]/(-64)*AscaleV/1000;
					}
				}
				*/
				// SPIKE DETECTION
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					//CHANNEL OUT OF LINEAR REGIME
					if (((vm[i][t+df]+4)%4096)<10) {
						if (A[i]<artT) {//reset only when it starts leaving the linear regime
							/*
							if (ACF) {
								for (int ii=0; ii<Math.Min(artT-A[i],6); ii++) {//is only for one step in the past...ignoring others
									SIprod [i, 12 - ii] -= (Aglobal[tA+dfSign]-Aglobal[tA+2*dfSign]) * (vm [i] [t + (6 - ii) * df] - vm [i] [t + (3 - ii) * df]) / Ascale;
								}
							}
							*/
							for (int ij=0; ij<dtMx; ij++) {
								Qdiff [i, ij] = 0;
							}
							A[i]=artT;
						}
					}
					//DEFAULT OPERATIONS
					else if (A[i]==0) {
						QmPreD [i] += (Avgs3[i] - QmPre [i]-QmPreD [i]) /(TauFilt);
						/*
						if (ACF) {
							SqIv [i] += (vm [i] [t + df] - vm [i] [t - 2 * df]) * (vm [i] [t + df] - vm [i] [t - 2 * df]);
							for (int iii=-6; iii<7; iii++) {
								SIprod [i, iii + 6] += Aglobaldiff[tA] * (vm [i] [t + (1 + iii) * df] - vm [i] [t - (2 - iii) * df]) / Ascale;
								//t-8...t+7
							}
						}
						*/
						Vbias[i]=FVbias[i]*AglobalSdiff[tA]/Sampling;
						//vm[i][t-df]+vm[i][t]+vm[i][t+df]
						Qdiff [i, dt]= (QmPre [i]+QmPreD [i])-Qm[i]-Vbias[i];//difference between ADC counts and Qm
						if ((AglobalSdiff[tA] * (Qdiff [i, dt] - Qdiff [i, dtPre])) != 0) {
							if ((AglobalSdiff[tA] > 0) == ((Qdiff [i, dt] - Qdiff [i, dtPre]) > 0)) {//(SIp[i]>0) {
								FVbias [i]++;
							} else {
								FVbias [i]--;
							}//Qdiff negative-->Ascale!!!;
						}
						//FVsbias[i]+=FVbias[i];//want to see whether this matches the correlation structure
						//Vsqbias[i]+=FVbias[i]*FVbias[i]/1000;
						//UPDATE Qm and Qd
						if (Qdiff[i,dt]>0) {
							if (Qdiff[i,dt]>Qd[i]) {
								Qm[i]+=Qd[i]/Tau_m0;
								if  (Qdiff[i,dt]<(5*Qd[i])) {
									Qd[i]++;
								}
								else if ((Qd[i]>Qdmin) & (Qdiff[i,dt]>(6*Qd[i]))) {
									Qd[i]--;
								}
							}
							else if (Qd[i]>Qdmin){//set a minimum level for Qd
								Qd[i]--;
							}
						}
						else if (Qdiff[i,dt]<-Qd[i]){
							Qm[i]-=Qd[i]/Tau_m0/2;
						}
					}
					//AFTER CHANNEL WAS OUT OF LINEAR REGIME
					else {
						Qm[i]=(2*Qm[i]+(QmPre [i]+QmPreD [i])+2*Qd[i])/3;//update Qm  vm[i][t-df]+vm[i][t]+vm[i][t+df]
						if (A [i] == artT) {
							QmPreD[i] = 0;
						}
						A[i]--;
					}
					//do above
					if (Qdiff [i, dt] * AmpScale > Qmax [i, dtPre]) {
						Qmax [i, dt] = Qdiff [i, dt] * AmpScale;//overflow issues?
						QmaxE [i] = dtE;
					} else if (dtE == QmaxE [i]) {
						if (Qdiff [i, dt] > Qdiff [i, dtPre]) {
							Qmax [i, dt] = Qdiff [i, dt] * AmpScale;//overflow issues?
							//QmaxE [i] = dtE;//redundant
						} else {
							Qmax [i, dt] = Qdiff [i, dtPre] * AmpScale;//overflow issues?
							QmaxE [i] = dtEx;
						}
					} else {
						Qmax [i, dt] = Qmax [i, dtPre];
						/*
						//thresholding
						if ((QmaxE [i] / Qd [i]) > (threshold - 2)) {
							for (int ii=0; ii<4; ii++) {
								if (ChInd5c [i] [ii] > -1) {
									T4 [ChInd5c [i] [ii]] = true;
								}
								if ((QmaxE [i] / Qd [i]) > (threshold - 2)) {
								}
							}
						}
						*/
					}
					//Qmax[i,dt]=Math.Max((Qdiff[i,(dt+1)%2],Qdiff[i,dt]);
				}
				foreach (int i in ChInd4List) {//loop across channels
					if (Sl4[i]>0) {//Sl frames after peak value
						ChS4Min = 0;
						a4=0;
						for (int ii=1; ii<4; ii++) {
							if (A[ChInd4a[i][ii]]==0){
								if (Qmax [ChInd4a[i][ChS4Min],dt]/Qd [ChInd4a[i][ChS4Min]] > Qmax [ChInd4a[i][ii],dt]/Qd [ChInd4a[i][ii]]){
								//if (Qmax [ChInd4a[i][ChS4Min],dt] > Qmax [ChInd4a[i][ii],dt]) {
									a4 += Qmax [ChInd4a[i][ChS4Min],dt]/ Qd [ChInd4a[i][ChS4Min]];
									ChS4Min = ii;
								}
								else {
									a4 += Qmax [ChInd4a[i][ii],dt]/ Qd [ChInd4a[i][ii]];
								}
							}
						}
						a4 /= 3;
						//default
						Sl4[i]=(Sl4[i]+1)%(Slmax+1);// increment Sl[i]
						//check whether it repolarizes
						if (a4>AHPthr) {
							AHP4[i]=true;
						}
						if (Sl4[i]==(Slmax-Sln0)){
							for (int ii=0; ii<4;ii++) {
								if (ChInd5a[i][ii]>-1){
									if (ChInd5a[i][ii]<ChInd4a[i][0]){//have updated Sl4 already
										if (Sl4[ChInd5a[i][ii]]>(Slmax-2*Sln0)){
											if (Amp4[ChInd5a[i][ii]]<Amp4[i]){
												Sl4x[ChInd5a[i][ii]]=true;
											}
										}
									}
									else {
										if ((Sl4[ChInd5a[i][ii]]>(Slmax-2*Sln0-1)) & (Sl4[ChInd5a[i][ii]]<(Slmax-1))){
											if (Amp4[ChInd5a[i][ii]]<Amp4[i]){
												Sl4x[ChInd5a[i][ii]]=true;
											}
										}
									}
								}
							}
							for (int ii=0; ii<4;ii++) {
								if (ChInd4c[i][ii]>-1){
									if ((Sl5[ChInd4c[i][ii]]>(Slmax-2*Sln0-1)) & (Sl5[ChInd4c[i][ii]]<(Slmax-1))){
										if (Amp5[ChInd4c[i][ii]]<Amp4[i]){
											Sl5x[ChInd4c[i][ii]]=true;
										}
									}
								}
							}
						}
						//accept spikes after Slmax frames if...
						if ((Sl4[i]==Slmax) & (!Sl4x[i])) {// & (AHP4[i]<(Slmax-Slmin))
							if (AHP4[i]) {
								wX.WriteLine("{0} {1} {2} {3} {4}", i, (t0+t)/dfAbs-dfSign*(Slmax-1), Amp4[i], 1, Acal-Slmax+1);
								NSpikes4++;//to count spikes
								foreach (int ch in ChInd4a[i]){
									for (int jj=t/dfAbs-CutOffset; jj<tCut+t/dfAbs-CutOffset; jj++) {
										tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
									}
									//compute baseline
									for (int kk=0; kk<tQm0; kk++){
										tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
									}
									for (int kk=tQm0; kk<tQmf; kk++){
										tQmA[kk]=tShape[tQmX[kk]];
									}
									Array.Sort (tQmA);
									b = tQmA [tQmm]+tQmA [tQmm+1];
									//compute max. amplitude above baseline
									Lmx = Lsw*b;
									for (int kk=Li; kk<Lf;kk++){
										Lm = Lw [0] * tShape [kk];
										for (int kkk=1; kkk<Lspike; kkk++) {
											Lm += Lw [kkk] * tShape [kkk+kk];
										}
										if (Lm<Lmx){
											Lmx=Lm;
										}
									}
									wShapesX.Write("{0} {1} {2} {3} ", ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
									for (int jj=0; jj<tCut; jj++){
										wShapesX.Write("{0} ", tShape[jj]);
									}
								}
								foreach (int ch in ChInd4b[i]){
									if (ch>-1){
										for (int jj=t/dfAbs-CutOffset; jj<tCut+t/dfAbs-CutOffset; jj++) {
											tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
										}
										//compute baseline
										for (int kk=0; kk<tQm0; kk++){
											tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
										}
										for (int kk=tQm0; kk<tQmf; kk++){
											tQmA[kk]=tShape[tQmX[kk]];
										}
										Array.Sort (tQmA);
										b = tQmA [tQmm]+tQmA [tQmm+1];
										//compute max. amplitude above baseline
										Lmx = Lsw*b;
										for (int kk=Li; kk<Lf;kk++){
											Lm = Lw [0] * tShape [kk];
											for (int kkk=1; kkk<Lspike; kkk++) {
												Lm += Lw [kkk] * tShape [kkk+kk];
											}
											if (Lm<Lmx){
												Lmx=Lm;
											}
										}
										wShapesX.Write("{0} {1} {2} {3} ",ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
										/*
										for (int jj=0; jj<tCut; jj++){
											wShapesX.Write("{0} ",  tShape[jj]);
										}
										*/
									}
									else {
										wShapesX.Write("{0} {1} {2} {3} ", ch,0,0,0);
										/*
										for (int jj=0; jj<tCut; jj++){
											wShapesX.Write("0 ");
										}
										*/
									}
								}
							}
							else {
								wX.WriteLine("{0} {1} {2} {3} {4}", i, (t0+t)/dfAbs-dfSign*(Slmax-1), Amp4[i], 0, Acal-Slmax+1);
								NSpikes4++;//to count spikes
								foreach (int ch in ChInd4a[i]){
									for (int jj=t/dfAbs-CutOffset; jj<tCutLong+t/dfAbs-CutOffset; jj++) {
										tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
									}
									//compute baseline
									for (int kk=0; kk<tQm0; kk++){
										tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
									}
									for (int kk=tQm0; kk<tQmf; kk++){
										tQmA[kk]=tShape[tQmXLong[kk]];
									}
									Array.Sort (tQmA);
									b = tQmA [tQmm]+tQmA [tQmm+1];
									//compute max. amplitude above baseline
									Lmx = Lsw*b;
									for (int kk=Li; kk<Lf;kk++){
										Lm = Lw [0] * tShape [kk];
										for (int kkk=1; kkk<Lspike; kkk++) {
											Lm += Lw [kkk] * tShape [kkk+kk];
										}
										if (Lm<Lmx){
											Lmx=Lm;
										}
									}
									wShapesX.Write("{0} {1} {2} {3} ", ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
									for (int jj=0; jj<tCutLong; jj++){
										wShapesX.Write("{0} ", tShape[jj]);
									}
								}
								foreach (int ch in ChInd4b[i]){
									if (ch>-1){
										for (int jj=t/dfAbs-CutOffset; jj<tCutLong+t/dfAbs-CutOffset; jj++) {
											tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
										}
										//compute baseline
										for (int kk=0; kk<tQm0; kk++){
											tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
										}
										for (int kk=tQm0; kk<tQmf; kk++){
											tQmA[kk]=tShape[tQmXLong[kk]];
										}
										Array.Sort (tQmA);
										b = tQmA [tQmm]+tQmA [tQmm+1];
										//compute max. amplitude above baseline
										Lmx = Lsw*b;
										for (int kk=Li; kk<Lf;kk++){
											Lm = Lw [0] * tShape [kk];
											for (int kkk=1; kkk<Lspike; kkk++) {
												Lm += Lw [kkk] * tShape [kkk+kk];
											}
											if (Lm<Lmx){
												Lmx=Lm;
											}
										}
										wShapesX.Write("{0} {1} {2} {3} ",ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
										/*
										for (int jj=0; jj<tCutLong; jj++){
											wShapesX.Write("{0} ", tShape[jj]);
										}
										*/
									}
									else {
										wShapesX.Write("{0} {1} {2} {3} ", ch,0,0,0);
										/*
										for (int jj=0; jj<tCutLong; jj++){
											wShapesX.Write("0 ");
										}
										*/
									}
								}
							}
							wShapesX.WriteLine();
							Sl4[i]=0;
						}
						//check whether current ADC count is higher
						else if (Amp4[i]<a4) {
							Sl4[i]=1;//reset peak value
							Sl4x[i]=false;
							Amp4[i]=a4;
							AHP4[i]=false;//reset AHP
						}
					}
					//}
					//else if (ChInd4a[i][0]>-1){//should rather make a list for that
					//want one step to decide whichelectrodes to look at; have to do for whole array on single electrodes.
					//if (Math.Max(Math.Min(Qmax[ChInd4a[i][0],dt],Qmax[ChInd4a[i][3],dt]),Math.Min(Qmax[ChInd4a[i][1],dt],Qmax[ChInd4a[i][2],dt]))>(2000*AmpScale)){
					else if (Z4next[i]==0) {
						ChS4Min = 0;
						a4 = 0;
						for (int ii=1; ii<4; ii++) {
							if (A [ChInd4a [i] [ii]] == 0) {
								if (Qmax [ChInd4a[i][ChS4Min],dt] / Qd [ChInd4a [i] [ChS4Min]] > Qmax [ChInd4a[i][ii],dt] / Qd [ChInd4a [i] [ii]]) {
									//if (Qmax [ChInd4a[i][ChS4Min],dt] > Qmax [ChInd4a[i][ii],dt]) {
									a4 += Qmax [ChInd4a [i] [ChS4Min],dt] / Qd [ChInd4a [i] [ChS4Min]];
									ChS4Min = ii;
								} else {
									a4 += Qmax [ChInd4a [i] [ii],dt] / Qd [ChInd4a [i] [ii]];
								}
							}
						}
						a4 /= 3;
						//check for threshold crossings
						if (a4 > threshold) {
							Sl4 [i] = 1;
							Sl4x [i] = false;
							Amp4 [i] = a4;
							AHP4 [i] = false;
						} else if (a4 < threshold-2*AmpScale) {
							Z4next[i] = dtTMx;
						}
					} else if (Z4next[i]==1) {//check for previous one as well
						Z4next[i] --;
						ChS4Min = 0;
						a4 = 0;
						for (int ii=1; ii<4; ii++) {
							if (A [ChInd4a [i] [ii]] == 0) {
								if (Qmax [ChInd4a[i][ChS4Min],dt] / Qd [ChInd4a [i] [ChS4Min]] > Qmax [ChInd4a[i][ii],dt] / Qd [ChInd4a [i] [ii]]) {
									//if (Qmax [ChInd4a[i][ChS4Min],dt] > Qmax [ChInd4a[i][ii],dt]) {
									a4 += Qmax [ChInd4a [i] [ChS4Min],dt] / Qd [ChInd4a [i] [ChS4Min]];
									ChS4Min = ii;
								} else {
									a4 += Qmax [ChInd4a [i] [ii],dt] / Qd [ChInd4a [i] [ii]];
								}
							}
						}
						a4 /= 3;
						//check for previous threshold crossing
						if (a4 > threshold) {
							Sl4 [i] = 1;
							Sl4x [i] = false;
							Amp4 [i] = a4;
							AHP4 [i] = false;
							ChS4Min = 0;
							a4 = 0;
							for (int ii=1; ii<4; ii++) {
								if (A [ChInd4a [i] [ii]] == 0) {
									if (Qmax [ChInd4a [i] [ChS4Min], dtPre] / Qd [ChInd4a [i] [ChS4Min]] > Qmax [ChInd4a [i] [ii], dtPre] / Qd [ChInd4a [i] [ii]]) {
										//if (Qmax [ChInd4a[i][ChS4Min],dtPre] > Qmax [ChInd4a[i][ii],dtPre]) {
										a4 += Qmax [ChInd4a [i] [ChS4Min], dtPre] / Qd [ChInd4a [i] [ChS4Min]];
										ChS4Min = ii;
									} else {
										a4 += Qmax [ChInd4a [i] [ii], dtPre] / Qd [ChInd4a [i] [ii]];
									}
								}
							}
							a4 /= 3;
							//check whether previous ADC count is higher
							if (Amp4 [i] < a4) {
								Sl4 [i] = 2;//reset peak value;
								Amp4 [i] = a4;
							}
						} else if (a4 < threshold -2*AmpScale) {
							Z4next [i] = dtTMx;
						}
						//check for previous threshold crossing (not sure whether necessary here, but wouldn't happen often)
						else {
							ChS4Min = 0;
							a4 = 0;
							for (int ii=1; ii<4; ii++) {
								if (A [ChInd4a [i] [ii]] == 0) {
									if (Qmax [ChInd4a [i] [ChS4Min], dtPre] / Qd [ChInd4a [i] [ChS4Min]] > Qmax [ChInd4a [i] [ii], dtPre] / Qd [ChInd4a [i] [ii]]) {
										//if (Qmax [ChInd4a[i][ChS4Min],dtPre] > Qmax [ChInd4a[i][ii],dtPre]) {
										a4 += Qmax [ChInd4a [i] [ChS4Min], dtPre] / Qd [ChInd4a [i] [ChS4Min]];
										ChS4Min = ii;
									} else {
										a4 += Qmax [ChInd4a [i] [ii], dtPre] / Qd [ChInd4a [i] [ii]];
									}
								}
							}
							a4 /= 3;
							if (a4 > threshold) {
								Sl4 [i] = 2;
								Sl4x [i] = false;
								Amp4 [i] = a4;
								AHP4 [i] = false;
							}
						}
					} else {
						Z4next[i] --;
					}
					//}
				}
				foreach (int i in ChInd5List) {//loop across channels
					if (Sl5[i]>0) {//Sl frames after peak value
						ChS5Min = Qmax [ChInd5a[i][0],dt]/Qd [ChInd5a[i][0]];
						a5=0;
						for (int ii=1; ii<4; ii++) {
							if (ChS5Min < Qmax [ChInd5a[i][ii],dt]/Qd [ChInd5a[i][ii]]) {
								a5 += Qmax [ChInd5a[i][ii],dt] / Qd [ChInd5a[i][ii]];
							}
							else {
								a5 += ChS5Min;
								ChS5Min=Qmax [ChInd5a[i][ii],dt] / Qd [ChInd5a[i][ii]];
							}
						}
						a5+=4*Qmax[ChInd5[i],dt]/ Qd [ChInd5[i]];
						a5 /= 7;
						//TREATMENT OF THRESHOLD CROSSINGS
						//default
						Sl5[i]=(Sl5[i]+1)%(Slmax+1);// increment Sl[i]
						//check whether it doesn't repolarize
						if (a5<AHPthr) {
							AHP5[i]=true;
						}
						if ((Sl5[i]==(Slmax-Sln0))){
							for (int ii=0; ii<4;ii++) {
								if (ChInd5a[i][ii]<ChInd5[i]){//have updated Sl5 already
									if (Sl5[ChInd5a[i][ii]]>(Slmax-2*Sln0)){
										if (Amp5[ChInd5a[i][ii]]<Amp5[i]){
											Sl5x[ChInd5a[i][ii]]=true;
										}
									}
								}
								else {
									if ((Sl5[ChInd5a[i][ii]]>(Slmax-2*Sln0-1)) & (Sl5[ChInd5a[i][ii]]<(Slmax-1))){
										if (Amp5[ChInd5a[i][ii]]<Amp5[i]){
											Sl5x[ChInd5a[i][ii]]=true;
										}
									}
								}
							}
							for (int ii=0; ii<4;ii++) {
								if (ChInd5c[i][ii]>-1) {
									if (Sl4[ChInd5c[i][ii]]>(Slmax-2*Sln0)){//have updated Sl4 already
										if (Amp4[ChInd5c[i][ii]]<Amp5[i]){
											Sl4x[ChInd5c[i][ii]]=true;
										}
									}
								}
							}
						}
						//accept spikes after Slmax frames if...
						if ((Sl5[i]==Slmax) & (!Sl5x[i])) {
							if (AHP5[i]) {
								w.WriteLine("{0} {1} {2} {3} {4}", i, (t0+t)/dfAbs-dfSign*(Slmax-1), Amp5[i], 1, Acal-Slmax+1);
								NSpikes5++;//to count spikes
								for (int jj=t/dfAbs-CutOffset; jj<tCut+t/dfAbs-CutOffset; jj++) {
									tShape[jj+CutOffset-t/dfAbs] = vm[ChInd5[i]][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ChInd5[i]]*Aglobaldiff[jj]/Sampling;
								}
								//compute baseline
								for (int kk=0; kk<tQm0; kk++){
									tQmA[kk]=(QmPreD[ChInd5[i]]-Qm[ChInd5[i]])*2/Ascale;
								}
								for (int kk=tQm0; kk<tQmf; kk++){
									tQmA[kk]=tShape[tQmX[kk]];
								}
								Array.Sort (tQmA);
								b = tQmA [tQmm]+tQmA [tQmm+1];
								//compute max. amplitude above baseline
								Lmx = Lsw*b;
								for (int kk=Li; kk<Lf;kk++){
									Lm = Lw [0] * tShape [kk];
									for (int kkk=1; kkk<Lspike; kkk++) {
										Lm += Lw [kkk] * tShape [kkk+kk];
									}
									if (Lm<Lmx){
										Lmx=Lm;
									}
								}
								wShapes.Write("{0} {1} {2} {3} ", ChInd5[i], ((QmPreD[ChInd5[i]]-Qm[ChInd5[i]])*2/Ascale), b, Lmx);
								for (int jj=0; jj<tCut; jj++){
									wShapes.Write("{0} ", tShape[jj]);
								}
								foreach (int ch in ChInd5a[i]){
									for (int jj=t/dfAbs-CutOffset; jj<tCut+t/dfAbs-CutOffset; jj++) {
										tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
									}
									//compute baseline
									for (int kk=0; kk<tQm0; kk++){
										tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
									}
									for (int kk=tQm0; kk<tQmf; kk++){
										tQmA[kk]=tShape[tQmX[kk]];
									}
									Array.Sort (tQmA);
									b = tQmA [tQmm]+tQmA [tQmm+1];
									//compute max. amplitude above baseline
									Lmx = Lsw*b;
									for (int kk=Li; kk<Lf;kk++){
										Lm = Lw [0] * tShape [kk];
										for (int kkk=1; kkk<Lspike; kkk++) {
											Lm += Lw [kkk] * tShape [kkk+kk];
										}
										if (Lm<Lmx){
											Lmx=Lm;
										}
									}
									wShapes.Write("{0} {1} {2} {3} ", ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
									for (int jj=0; jj<tCut; jj++){
										wShapes.Write("{0} ", tShape[jj]);
									}
								}
								foreach (int ch in ChInd5b[i]){
									if (ch>-1){
										for (int jj=t/dfAbs-CutOffset; jj<tCut+t/dfAbs-CutOffset; jj++) {
											tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
										}
										//compute baseline
										for (int kk=0; kk<tQm0; kk++){
											tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
										}
										for (int kk=tQm0; kk<tQmf; kk++){
											tQmA[kk]=tShape[tQmX[kk]];
										}
										Array.Sort (tQmA);
										b = tQmA [tQmm]+tQmA [tQmm+1];
										//compute max. amplitude above baseline
										Lmx = Lsw*b;
										for (int kk=Li; kk<Lf;kk++){
											Lm = Lw [0] * tShape [kk];
											for (int kkk=1; kkk<Lspike; kkk++) {
												Lm += Lw [kkk] * tShape [kkk+kk];
											}
											if (Lm<Lmx){
												Lmx=Lm;
											}
										}
										wShapes.Write("{0} {1} {2} {3} ", ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
										/*
										for (int jj=0; jj<tCut; jj++){
											wShapes.Write("{0} ", tShape[jj]);
										}
										*/
									}
									else {
										wShapes.Write("{0} {1} {2} {3} ", ch,0,0,0);
										/*
										for (int jj=0; jj<tCut; jj++){
											wShapes.Write("0 ");
										}
										*/
									}
								}
							}
							else {
								w.WriteLine("{0} {1} {2} {3} {4}", i, (t0+t)/dfAbs-dfSign*(Slmax-1), Amp5[i], 0, Acal-Slmax+1);
								NSpikes5++;//to count spikes
								for (int jj=t/dfAbs-CutOffset; jj<tCutLong+t/dfAbs-CutOffset; jj++) {
									tShape[jj+CutOffset-t/dfAbs] = vm[ChInd5[i]][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ChInd5[i]]*Aglobaldiff[jj]/Sampling;
								}
								//compute baseline
								for (int kk=0; kk<tQm0; kk++){
									tQmA[kk]=(QmPreD[ChInd5[i]]-Qm[ChInd5[i]])*2/Ascale;
								}
								for (int kk=tQm0; kk<tQmf; kk++){
									tQmA[kk]=tShape[tQmXLong[kk]];
								}
								Array.Sort (tQmA);
								b = tQmA [tQmm]+tQmA [tQmm+1];
								//compute max. amplitude above baseline
								Lmx = Lsw*b;
								for (int kk=Li; kk<Lf;kk++){
									Lm = Lw [0] * tShape [kk];
									for (int kkk=1; kkk<Lspike; kkk++) {
										Lm += Lw [kkk] * tShape [kkk+kk];
									}
									if (Lm<Lmx){
										Lmx=Lm;
									}
								}
								wShapes.Write("{0} {1} {2} {3} ", ChInd5[i], ((QmPreD[ChInd5[i]]-Qm[ChInd5[i]])*2/Ascale), b, Lmx);
								for (int jj=0; jj<tCutLong; jj++){
									wShapes.Write("{0} ", tShape[jj]);
								}
								foreach (int ch in ChInd5a[i]){
									for (int jj=t/dfAbs-CutOffset; jj<tCutLong+t/dfAbs-CutOffset; jj++) {
										tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
									}
									//compute baseline
									for (int kk=0; kk<tQm0; kk++){
										tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
									}
									for (int kk=tQm0; kk<tQmf; kk++){
										tQmA[kk]=tShape[tQmXLong[kk]];
									}
									Array.Sort (tQmA);
									b = tQmA [tQmm]+tQmA [tQmm+1];
									//compute max. amplitude above baseline
									Lmx = Lsw*b;
									for (int kk=Li; kk<Lf;kk++){
										Lm = Lw [0] * tShape [kk];
										for (int kkk=1; kkk<Lspike; kkk++) {
											Lm += Lw [kkk] * tShape [kkk+kk];
										}
										if (Lm<Lmx){
											Lmx=Lm;
										}
									}
									wShapes.Write("{0} {1} {2} {3} ", ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
									for (int jj=0; jj<tCutLong; jj++){
										wShapes.Write("{0} ", tShape[jj]);
									}
								}
								foreach (int ch in ChInd5b[i]){
									if (ch>-1){
										for (int jj=t/dfAbs-CutOffset; jj<tCutLong+t/dfAbs-CutOffset; jj++) {
											tShape[jj+CutOffset-t/dfAbs] = vm[ch][jj*dfAbs]*2 -Aglobal[jj] -FVbias[ch]*Aglobaldiff[jj]/Sampling;
										}
										//compute baseline
										for (int kk=0; kk<tQm0; kk++){
											tQmA[kk]=(QmPreD[ch]-Qm[ch])*2/Ascale;
										}
										for (int kk=tQm0; kk<tQmf; kk++){
											tQmA[kk]=tShape[tQmXLong[kk]];
										}
										Array.Sort (tQmA);
										b = tQmA [tQmm]+tQmA [tQmm+1];
										//compute max. amplitude above baseline
										Lmx = Lsw*b;
										for (int kk=Li; kk<Lf;kk++){
											Lm = Lw [0] * tShape [kk];
											for (int kkk=1; kkk<Lspike; kkk++) {
												Lm += Lw [kkk] * tShape [kkk+kk];
											}
											if (Lm<Lmx){
												Lmx=Lm;
											}
										}
										wShapes.Write("{0} {1} {2} {3} ", ch, ((QmPreD[ch]-Qm[ch])*2/Ascale), b, Lmx);
										/*
										for (int jj=0; jj<tCutLong; jj++){
											wShapes.Write("{0} ", tShape[jj]);
										}
										*/
									}
									else {
										wShapes.Write("{0} {1} {2} {3} ", ch,0,0,0);
										/*
										for (int jj=0; jj<tCutLong; jj++){
											wShapes.Write("0 ");
										}
										*/
									}
								}
							}
							wShapes.WriteLine();
							Sl5[i]=0;
						}
						//check whether current ADC count is higher
						else if (Amp5[i]<a5) {
							Sl5[i]=1;//reset peak value
							Sl5x[i]=false;
							Amp5[i]=a5;
							AHP5[i]=false;//reset AHP
						}
						//}
					}
					//else if (ChInd5[i]>-1){
					//if (Qmax[ChInd5[i],dt]>(thr5*Qd[i])){
					else if (Z5next [i] == 0) {
						ChS5Min = Qmax [ChInd5a [i] [0], dt] / Qd [ChInd5a [i] [0]];
						a5 = 0;
						for (int ii=1; ii<4; ii++) {
							if (ChS5Min < Qmax [ChInd5a [i] [ii], dt] / Qd [ChInd5a [i] [ii]]) {
								a5 += Qmax [ChInd5a [i] [ii], dt] / Qd [ChInd5a [i] [ii]];
							} else {
								a5 += ChS5Min;
								ChS5Min = Qmax [ChInd5a [i] [ii], dt] / Qd [ChInd5a [i] [ii]];
							}
						}
						a5 += 4 * Qmax [ChInd5 [i], dt] / Qd [ChInd5 [i]];
						a5 /= 7;
						//check for threshold crossings
						if (a5 > threshold) {
							Sl5 [i] = 1;
							Sl5x [i] = false;
							Amp5 [i] = a5;
							AHP5 [i] = false;
						} else if (a5 < threshold -2*AmpScale) {
							Z5next [i] = dtTMx;
						}
					} else if (Z5next [i] == 1) {//check for previous one as well
						Z5next [i] --;
						ChS5Min = Qmax [ChInd5a [i] [0], dt] / Qd [ChInd5a [i] [0]];
						a5 = 0;
						for (int ii=1; ii<4; ii++) {
							if (ChS5Min < Qmax [ChInd5a [i] [ii], dt] / Qd [ChInd5a [i] [ii]]) {
								a5 += Qmax [ChInd5a [i] [ii], dt] / Qd [ChInd5a [i] [ii]];
							} else {
								a5 += ChS5Min;
								ChS5Min = Qmax [ChInd5a [i] [ii], dt] / Qd [ChInd5a [i] [ii]];
							}
						}
						a5 += 4 * Qmax [ChInd5 [i], dt] / Qd [ChInd5 [i]];
						a5 /= 7;
						//check for previous threshold crossing
						if (a5 > threshold) {
							Sl5 [i] = 1;
							Sl5x [i] = false;
							Amp5 [i] = a5;
							AHP5 [i] = false;
							ChS5Min = Qmax [ChInd5a [i] [0], dtPre] / Qd [ChInd5a [i] [0]];
							a5 = 0;
							for (int ii=1; ii<4; ii++) {
								if (ChS5Min < Qmax [ChInd5a [i] [ii], dtPre] / Qd [ChInd5a [i] [ii]]) {
									a5 += Qmax [ChInd5a [i] [ii], dtPre] / Qd [ChInd5a [i] [ii]];
								} else {
									a5 += ChS5Min;
									ChS5Min = Qmax [ChInd5a [i] [ii], dtPre] / Qd [ChInd5a [i] [ii]];
								}
							}
							a5 += 4 * Qmax [ChInd5 [i], dtPre] / Qd [ChInd5 [i]];
							a5 /= 7;
							//check for previous threshold crossing
							if (Amp5 [i] < a5) {
								Sl5 [i] = 2;
								Amp5 [i] = a5;
							}
						} else if (a5 < threshold -2*AmpScale) {
							Z5next [i] = dtTMx;
						}
						//check for previous threshold crossing (not sure whether necessary here, but wouldn't happen often)
						else {
							ChS5Min = Qmax [ChInd5a [i] [0], dtPre] / Qd [ChInd5a [i] [0]];
							a5 = 0;
							for (int ii=1; ii<4; ii++) {
								if (ChS5Min < Qmax [ChInd5a [i] [ii], dtPre] / Qd [ChInd5a [i] [ii]]) {
									a5 += Qmax [ChInd5a [i] [ii], dtPre] / Qd [ChInd5a [i] [ii]];
								} else {
									a5 += ChS5Min;
									ChS5Min = Qmax [ChInd5a [i] [ii], dtPre] / Qd [ChInd5a [i] [ii]];
								}
							}
							a5 += 4 * Qmax [ChInd5 [i], dtPre] / Qd [ChInd5 [i]];
							a5 /= 7;
							//check for previous threshold crossing
							if (a5 > threshold) {
								Sl5 [i] = 2;
								Sl5x [i] = false;
								Amp5 [i] = a5;
								AHP5 [i] = false;
							}
						}
					} else {
						Z5next[i] --;
					}
					//}
				}
			}
		}

		public void FinishDetection (short[][] vm, int skipLast)
		{//finish baseline estimate, write spikes in interval after last recalibration; close file
			if (df > 0) {
				for (int t=tf; t<df*(tInc-1)-skipLast; t+=df) {//loop over data, will be removed for an online algorithm
					for (int i=1; i<NChannels; i++) {//loop across channels
						Slice [i] = (vm [i] [t]) % 4095 + (vm [i] [t + df]) % 4095;
					}
					Array.Sort (Slice);
					wMean.WriteLine ("{0}", Slice [NChannels / 2]);
				}
			} else {
				for (int t=tf; t>0; t+=df) {//loop over data, will be removed for an online algorithm
					for (int i=1; i<NChannels; i++) {//loop across channels
						Slice [i] = (vm [i] [t]) % 4095 + (vm [i] [t + df]) % 4095;
					}
					Array.Sort (Slice);
					wMean.WriteLine ("{0}", Slice [NChannels / 2]);
				}
			}
			wMean.WriteLine ("{0}", Slice [NChannels / 2]);
			wInfo.WriteLine ("#Sum(squared global fluctuations):");
			wInfo.WriteLine ("{0}", SqIglobal);
			if (ACF) {
				wInfo.WriteLine ("#Sum(squared channel fluctuations):");
				for (int i=0; i<NChannels; i++) {//loop across channels
					wInfo.Write ("{0} ", SqIv[i]);
				}
				wInfo.WriteLine();
				wInfo.WriteLine ("#Sum(product of channel and global fluctuations):");
				for (int ii=0; ii<13; ii++) {//loop across timelags
					for (int i=0; i<NChannels; i++) {//loop across channels
						wInfo.Write ("{0} ", SIprod [i, ii]);
					}
					wInfo.WriteLine ();
				}
			}
			wInfo.WriteLine ("#Sum(avg. deviations from global fluctuations):");
			for (int i=0; i<NChannels; i++) {//loop across channels
				wInfo.Write ("{0} ", FVsbias[i]);
			}
			wInfo.WriteLine();
			wInfo.WriteLine ("#Sum(avg. squared deviations from global fluctuations):");
			for (int i=0; i<NChannels; i++) {//loop across channels
				wInfo.Write ("{0} ", Vsqbias[i]);
			}
			wInfo.WriteLine();
			wInfo.WriteLine("#Number of spikes (4 channel):");
			wInfo.WriteLine("{0} ",NSpikes4);//to count spikes
			wInfo.WriteLine("#Number of spikes (5 channel):");
			wInfo.WriteLine("{0} ",NSpikes5);//to count spikes
			w.Flush();
			w.Close();
			fs.Close();
			wShapes.Flush();
			wShapes.Close();
			fsShapes.Close();
			wX.Flush();
			wX.Close();
			fsX.Close();
			wShapesX.Flush();
			wShapesX.Close();
			fsShapesX.Close();
			wInfo.Flush();
			wInfo.Close();
			fsInfo.Close();
			wMean.Flush();
			wMean.Close();
			fsMean.Close();
		}
	}
	public class MainClass{
		[STAThread]
		static public void Main ()
		{
			// Initialize an instance of the BrwRdr class
			BrwRdr brwRdr = new BrwRdr ();
			// you must use a valid full file path
			
			// Create and display a fileChooserDialog
			OpenFileDialog openFileDialog1 = new OpenFileDialog ();
			//openFileDialog1.InitialDirectory = "c:\\";
			openFileDialog1.Filter = "brw files (*.brw)|*.brw|All files (*.*)|*.*";
			openFileDialog1.FilterIndex = 1;
			openFileDialog1.RestoreDirectory = true;
			openFileDialog1.InitialDirectory = Directory.GetCurrentDirectory ();

			if (openFileDialog1.ShowDialog () != DialogResult.OK) {
				throw new Exception ("File dialog error");
			}
			// Open the file for reading.
			brwRdr.Open (openFileDialog1.FileName);
			// the number of frames of the recording
			long nFrames = brwRdr.RecNFrames;
			// the duration in seconds of the recording
			double nSec = brwRdr.RecDuration / 1000;
			// the sampling frequency
			double sfd = (double)brwRdr.SamplingRate;
			int sf = (int)sfd;
			Console.WriteLine ("# Number of frames: {0}", nFrames);
			Console.WriteLine ("# Duration (s): {0}", nSec);
			Console.WriteLine ("# Sampling rate: {0}", sfd);

			//list of channels
			ChCoord[] Channels = brwRdr.GetRecChs (StreamType.Raw);
			Console.WriteLine ("# Number of recorded channels: {0}", Channels.Length);
			int[] Indices = new int[Channels.Length];
			for (int i=0; i<Channels.Length; i++) {
				Indices [i] = (Channels [i].Col - 1) + 64 * (Channels [i].Row - 1);
				Console.WriteLine ("{0}", Indices [i]);
			}
			// how many frames to analyze
			//long nDumpFrames = nFrames;
			Detection SpkD;
			SpkD = new Detection ();

			//set initial parameter
			int[] dfTI = SpkD.SetInitialParams (nFrames, nSec, sf, sfd, Channels.Length, Indices);

			SaveFileDialog saveFileDialog1 = new SaveFileDialog ();
			saveFileDialog1.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
			saveFileDialog1.FilterIndex = 1;
			saveFileDialog1.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Spikes";
			saveFileDialog1.DefaultExt = "txt";
			saveFileDialog1.InitialDirectory = Directory.GetCurrentDirectory ();
			saveFileDialog1.Title = "Save Spikes As";
				
			if (saveFileDialog1.ShowDialog () != DialogResult.OK) {
				throw new Exception ("File dialog error");
			}
			SaveFileDialog saveFileDialog2 = new SaveFileDialog ();
			saveFileDialog2.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
			saveFileDialog2.FilterIndex = 1;
			saveFileDialog2.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Shapes";
			saveFileDialog2.DefaultExt = "txt";
			saveFileDialog2.Title = "Save Shapes As";

			if (saveFileDialog2.ShowDialog () != DialogResult.OK) {
				throw new Exception ("File dialog error");
			}
			SaveFileDialog saveFileDialog5 = new SaveFileDialog ();
			saveFileDialog5.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
			saveFileDialog5.FilterIndex = 1;
			saveFileDialog5.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_SpikesX";
			saveFileDialog5.DefaultExt = "txt";
			saveFileDialog5.Title = "Save Interpolated Shapes As";

			if (saveFileDialog5.ShowDialog () != DialogResult.OK) {
				throw new Exception ("File dialog error");
			}
			SaveFileDialog saveFileDialog6 = new SaveFileDialog ();
			saveFileDialog6.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
			saveFileDialog6.FilterIndex = 1;
			saveFileDialog6.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_ShapesX";
			saveFileDialog6.DefaultExt = "txt";
			saveFileDialog6.Title = "Save Interpolated Shapes As";

			if (saveFileDialog6.ShowDialog () != DialogResult.OK) {
				throw new Exception ("File dialog error");
			}
			SaveFileDialog saveFileDialog3 = new SaveFileDialog ();
			saveFileDialog3.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
			saveFileDialog3.FilterIndex = 1;
			saveFileDialog3.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Info";
			saveFileDialog3.DefaultExt = "txt";
			saveFileDialog3.Title = "Save Further Information As";

			if (saveFileDialog3.ShowDialog () != DialogResult.OK) {
				throw new Exception ("File dialog error");
			}
			SaveFileDialog saveFileDialog4 = new SaveFileDialog ();
			saveFileDialog4.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
			saveFileDialog4.FilterIndex = 1;
			saveFileDialog4.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Avg";
			saveFileDialog4.DefaultExt = "txt";
			saveFileDialog4.Title = "Save Average Voltage Increments As";

			if (saveFileDialog4.ShowDialog () != DialogResult.OK) {
				throw new Exception ("File dialog error");
			}

			// open the output file
			SpkD.openSpikeFile (saveFileDialog1.FileName);
			saveFileDialog1.Dispose ();
			SpkD.openShapeFile (saveFileDialog2.FileName);
			saveFileDialog2.Dispose ();
			SpkD.openSpikeXFile (saveFileDialog5.FileName);
			saveFileDialog5.Dispose ();
			SpkD.openShapeXFile (saveFileDialog6.FileName);
			saveFileDialog6.Dispose ();
			SpkD.openInfoFile (saveFileDialog3.FileName);
			saveFileDialog3.Dispose ();
			SpkD.openMeanFile (saveFileDialog4.FileName);
			saveFileDialog4.Dispose ();

			// measure execution time
			var sw = new Stopwatch ();
			sw.Start ();
			//const int NChannels = 4096;
			//const int df = 2;//has to be changed in Detection class as well
			if (dfTI [0] > 0) {
				long t1 = 0;
				long tInc = dfTI [2];//has to be changed in Detection class as well
				//int tCut = (sf / 501 + sf / 1000) / dfTI[0] + 6;
				//int CutOffset = (sf / 1000) / df + 6;
				for (long t0=0; t0<Math.Min(200*(dfTI[1]),nFrames-tInc); t0+=dfTI[1]) {
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
					SpkD.InitialEstimation (vm, t0);
				}
				for (long t0=0; t0<dfTI[1]; t0+=dfTI[1]) {
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
					SpkD.StartDetection (vm, t0, nFrames, nSec, sfd, Indices);
					SpkD.Iterate (vm, t0);
					t1 += dfTI [1];
				}
				for (long t0=dfTI[1]; t0<nFrames-tInc; t0+=dfTI[1]) {
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
					SpkD.Iterate (vm, t0);
					t1 += dfTI [1];
				}
				if (t1 < nFrames - tInc + dfTI [1] - 1) {
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
					SpkD.skipLastReverse ((int)(tInc - nFrames + t1));
					SpkD.Iterate (vm, t1);
				}
				short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
				SpkD.FinishDetection (vmx, (int)(tInc - nFrames + t1));
			} else {
				long t1 = nFrames;
				long tInc = dfTI [2];
				//int tCut = -(sf / 501 + sf / 1000) / df + 6 +8;
				//int CutOffset = -(sf / 1000) / df + 6;
				for (long t0=nFrames; t0>Math.Max(tInc,nFrames-200*dfTI[1]); t0-= dfTI[1]) {
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
					SpkD.InitialEstimation (vm, t0 - tInc);
				}
				for (long t0=nFrames; t0>nFrames-dfTI[1]; t0-=dfTI[1]) {
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0-tInc, tInc);
					SpkD.StartDetection (vm, t0-tInc, nFrames, nSec, sfd, Indices);
					SpkD.Iterate (vm, t0-tInc);
					t1 -= dfTI [1];
				}
				for (long t0=nFrames-dfTI[1]; t0>tInc; t0-= dfTI[1]) {
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
					SpkD.Iterate (vm, t0 - tInc);
					t1 -= dfTI [1];
				}
				if (t1 > tInc - dfTI [1] + 1) {
					SpkD.skipLastReverse ((int)(tInc - t1));
					short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
					SpkD.Iterate (vm, 0);
				}
				short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
				SpkD.FinishDetection (vmx, (int)(tInc - t1));
			}
			sw.Stop ();
			Console.WriteLine ("Elapsed time: {0}", sw.Elapsed); // TimeSpan
			Console.WriteLine ("Milliseconds/frame: {0}", sw.Elapsed.TotalMilliseconds / nFrames); // TimeSpan
		}
	}
}





