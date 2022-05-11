#include <iostream>
#include <fstream>
#include <string>

#include <TFile.h>
//#include <TMath.h>
//#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>


// 1st Method
// hist is a pointer to https://root.cern.ch/doc/master/classTH1D.html
void get_peak_centres_method1(TH1F* hist, Double_t* peak_Q_values, Double_t* secondary_Q_values, Double_t* peak_centres, Double_t* secondary_centres) {

	std::string ver = "\nget_peak_centres_simple_method v1\n";
	std::cout << ver << std::endl;

	// Results array [0-3][0-1] 0=bin centre, 1=count
	Double_t max_peaks[4][2] = { {0,0},{0,0},{0,0},{0,0} };

	int count;
	int min;

	int min_index;
	// No of channels to lookahead to for searh of max
	int peak_search_bins = 6;
	// Ignore everything below the specified channel
	int start_channel = 280;

	// Loop through the strip data looking for max peaks
	int i = start_channel;
	while (i < 4096) {

		count = hist->GetBinContent(i);
		min = std::min({ max_peaks[0][1], max_peaks[1][1], max_peaks[2][1], max_peaks[3][1] });

		// Do we have a new max peak
		if (count > min) {

			// Replace the current lowest bin with the new max bin
			if (max_peaks[0][1] == min) {
				min_index = 0;
			}
			else if (max_peaks[1][1] == min) {
				min_index = 1;
			}
			else if (max_peaks[2][1] == min) {
				min_index = 2;
			}
			else {
				min_index = 3;
			}

			for (int j = 0; j < peak_search_bins; j++) {
				count = hist->GetBinContent(i + j);
				if (count > max_peaks[min_index][1]) {
					max_peaks[min_index][0] = hist->GetBinCenter(i + j);  // bin centre
					max_peaks[min_index][1] = count;                      // count
				}
			}

			// Having found a peak continue searching after the last compared bin
			i = i + peak_search_bins - 1;
		} // end if

		i++;
	} // end loop

	// index = 0 bin centre, 1 count
	i = 0;
	// Having found the max peaks return an array of channel bin centres
	peak_centres[0] = max_peaks[0][i];
	peak_centres[1] = max_peaks[1][i];
	peak_centres[2] = max_peaks[2][i];
	peak_centres[3] = max_peaks[3][i];

	return;
}


// 2nd Method
// hist is a pointer to https://root.cern.ch/doc/master/classTH1D.html
void get_peak_centres_method2(TH1F* hist, Double_t* peak_Q_values, Double_t* secondary_Q_values, Double_t* peak_centres, Double_t* secondary_centres) {

	std::string ver = "\nget_peak_centres_using_fit_gauss1 v1\n";
	std::cout << ver << std::endl;

	// Results array [0-3][0-3] 0=bin#, 1=count, 2=mean, 3=sigma
	Double_t max_peaks[4][4] = { {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0} };

	int count;
	int min;

	int min_index;
	// No of channels to read-ahead to for search of max
	int peak_search_bins = 8;
	// Ignore everything below the specified channel
	int start_channel = 280;

	// Loop through the strip data looking for max peaks
	int i = start_channel;
	while (i < 4096) {

		count = hist->GetBinContent(i);
		min = std::min({ max_peaks[0][1], max_peaks[1][1], max_peaks[2][1], max_peaks[3][1] });

		if (count > min) {

			if (max_peaks[0][1] == min) {
				min_index = 0;
			}
			else if (max_peaks[1][1] == min) {
				min_index = 1;
			}
			else if (max_peaks[2][1] == min) {
				min_index = 2;
			}
			else {
				min_index = 3;
			}

			for (int j = 0; j < peak_search_bins; j++) {
				count = hist->GetBinContent(i + j);
				if (count > max_peaks[min_index][1]) {
					max_peaks[min_index][0] = hist->GetBinCenter(i + j);  // bin centre
					max_peaks[min_index][1] = count;                      // count
				}
			}

			// Fit a Gaussian to the peak
			hist->Fit("gaus", "Q", "", max_peaks[min_index][0] - (peak_search_bins / 2), max_peaks[min_index][0] + (peak_search_bins / 2));
			TF1* gausfit = (TF1*)hist->GetFunction("gaus");
			// hist->Draw();

			// GetParameter(0 Constant, 1 Mean, 2 Sigma)
			// Gaussian Fit - mean
			max_peaks[min_index][2] = gausfit->GetParameter(1);
			// Gaussian Fit - sigma
			max_peaks[min_index][3] = gausfit->GetParameter(2);

			std::cout << "For bin centre = " + std::to_string(max_peaks[min_index][0]) << std::endl;
			std::cout << "\t- count = " + std::to_string(max_peaks[min_index][1]) << std::endl;
			std::cout << "\t- mean = " + std::to_string(max_peaks[min_index][2]) << std::endl;
			std::cout << "\t- sigma = " + std::to_string(max_peaks[min_index][3]) << std::endl;

			// Having found a peak continue searching from mean (peak) + 3*sigma (round up)
			i = std::round(max_peaks[min_index][2] + 3 * max_peaks[min_index][3] + 0.5);

		} // end if

		i++;
	} // end loop


	// index = 0 bin centre, 1 count, 2 mean, 3 sigma
	i = 0;
	// Having found the max peaks return an array of channel bin centres
	peak_centres[0] = max_peaks[0][i];
	peak_centres[1] = max_peaks[1][i];
	peak_centres[2] = max_peaks[2][i];
	peak_centres[3] = max_peaks[3][i];

	return;
}

// 3rd Method
// hist is a pointer to https://root.cern.ch/doc/master/classTH1D.html
void get_peak_centres_method3(TH1F* hist, Double_t* peak_Q_values, Double_t* secondary_Q_values, Double_t* peak_centres, Double_t* secondary_centres) {

	std::string ver = "\nget_peak_centres_using_fit_Q_values v1\n";
	std::cout << ver << std::endl;

	// Results array [0-3][0-3] 0=bin#, 1=count, 2=mean, 3=sigma
	Double_t max_peaks[4][4] = { {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0} };

	// No of channels to read-ahead to search of max
	int peak_search_bins = 20;
	// Ignore everything below the specified channel
	int start_channel = 280;

	double factor = 1 / 15.5;
	double offset = 230;
	
	std::cout << "Factor " << factor << "\t" << 1 / factor << std::endl;

	// Loop through the strip data looking for max peaks
	int i = start_channel;
	int p = 0;
	int bin;
	int count;
	int max;
	while (i < 4096) {

		// Start the search for each main peak based on its Q-value
		std::cout << "Got Q value " + std::to_string(peak_Q_values[p]) << std::endl;
		// Apply factor and offset to each Q-value
		i = peak_Q_values[p] * factor + offset;
		bin = 0;
		max = 0;

		std::cout << "Search starting at " + std::to_string(i) << std::endl;

		for (int j = 0; j < peak_search_bins; j++) {
			count = hist->GetBinContent(i + j);
			if (count > max) {
				bin = i + j;
				max = count;
			}
		}

		// Record bin centre & count
		max_peaks[p][0] = hist->GetBinCenter(bin);  // bin centre
		max_peaks[p][1] = count;                    // count

		// Having searched for a peak continue searching from end of search
		i = i + peak_search_bins;

		p++;

		// Stop search after we've found the 4th peak
		if (p == 4) {
			break;
		}

	} // end loop


	// index = 0 bin centre, 1 count
	i = 0;
	// Having found the max peaks return an array of channel bin centres
	peak_centres[0] = max_peaks[0][i];
	peak_centres[1] = max_peaks[1][i];
	peak_centres[2] = max_peaks[2][i];
	peak_centres[3] = max_peaks[3][i];

	return;
}

// hist is a pointer to https://root.cern.ch/doc/master/classTH1D.html
void get_peak_centres_method4(TH1F* hist, Double_t* peak_Q_values, Double_t* secondary_Q_values, Double_t* peak_centres, Double_t* secondary_centres) {

	std::string ver = "\nget_peak_centres_using_fit_first_then_offset v2\n";
	std::cout << ver << std::endl;

	double factor = 1 / 15.5;
//	std::cout << "Factor " << factor << "\t" << 1 / factor << std::endl;

	// Bin of 1st max peak 
	int max_bin;
	
	// Variables used to process each remaining peaks
	// No of channels to read-ahead when searching for max
	int peak_search_bins;
	// Used to track max bin
	int bin;
	// Used to track count
	int count;
	// Used to track max count
	int max;

	// Bin to start peak search
	int start;

	// Find the bin of 1st max peak 
//	max_bin = hist->GetBinCenter(hist->GetMaximumBin());

	// Loop through the strip data looking for the 'max' peak
	// Ignore everything below the specified channel
	int start_channel = 320;
	int end_channel = 450;
	int max_count = 0;
	int i = start_channel;
	while (i < 4096) {

		count = hist->GetBinContent(i);
		if (count > max_count) {
			max_count = count;
			max_bin = i;
		} // end if

		if (i==end_channel) {
			break;
		}
		
		i++;
	} // end loop

// OLD	max_bin = hist->GetMaximumBin();
//	max_bin = hist->GetMaximumBin();

	peak_centres[0] = max_bin;
	std::cout << "Peak 0 - max bin at centre " << max_bin << std::endl;

	// Loop through remaining main peak Q-values
	peak_search_bins = 12;
	for (int i = 1; i < 4; i++) {

		bin = count = max = 0;
		start = hist->GetBinCenter(max_bin + ((peak_Q_values[i] - peak_Q_values[0]) * factor)) - peak_search_bins / 2;

		for (int j = start; j < start + peak_search_bins; j++) {
			count = hist->GetBinContent(j);
			if (count > max) {
				bin = j;
				max = count;
			}
		}
//		peak_centres[i] = hist->GetBinCenter(bin);
		peak_centres[i] = bin;
		std::cout << "Peak " << i << " - started at " << start << " ended at " << start + peak_search_bins << " found max " << max << " at " << bin << " centre " << peak_centres[i] << std::endl;
	}

	// Loop through fine structure Q values
	peak_search_bins = 4;
	for (int i = 0; i < 4; i++) {

		// Use different read-ahead for peaks 2 & 3
		if (i == 2) {
			peak_search_bins = 6;
		}

		bin = count = max = 0;
		start = hist->GetBinCenter(max_bin + ((secondary_Q_values[i] - peak_Q_values[0]) * factor)) - peak_search_bins / 2;

		for (int j = start; j < start + peak_search_bins; j++) {
			count = hist->GetBinContent(j);
			if (count > max) {
				bin = j;
				max = count;
			}
		}
//		secondary_centres[i] = hist->GetBinCenter(bin);
		secondary_centres[i] = bin;
		std::cout << "Secondary " << i << " - started at " << start << " ended at " << start + peak_search_bins << " found max " << max << " at " << bin << " centre " << secondary_centres[i] << std::endl;
	}
	return;
}


// *************************
// END OF PEAK SERCH METHODS
// *************************


// Background function
double fitBackground(double* x, double* par) {
	// y = c + mx
	return par[0] + par[1] * x[0];
}

// Peak function - Gaussian
double fitPeak(double* x, double* par) {
	// 0 - height, 1 - centre, 2 - std dev (sigma)
//	return par[0] * exp(-0.5 * std::pow((x[0] - par[1]) / par[2], 2));
	double sigma[4] = {1, 1, 1, 1};
//	return par[0] * TMath::Gaus(x[0], par[1], 3);
	return par[0] * TMath::Gaus(x[0], par[1], par[2]);
}

// Peak function - Gaussian
double fitPeak2(double* x, double* par, int sigma_index) {
	// 0 - height, 1 - centre, 2 - std dev (sigma)
//	return par[0] * exp(-0.5 * std::pow((x[0] - par[1]) / 2, 2));
//	double sigma[4] = {1.2, 1.2, 1.4, 2};
	double sigma[4] = {0.8, 1.2, 1.2, 0.25};
//	double sigma[4] = {0.25, 0.25, 0.25, 0.25};
//	double sigma[4] = {0.491233, 0.474015, 0.473551, 0.483108};
	return par[0] * TMath::Gaus(x[0], par[1], sigma[sigma_index]);
}

// Sum of background and peak function
double fitFunction(double* x, double* par) {
	return fitBackground(x, &par[0]) \
		+ fitPeak(x, &par[2]) + fitPeak(x, &par[5]) + fitPeak(x, &par[8]) + fitPeak(x, &par[11]) \
		+ fitPeak(x, &par[14]) + fitPeak(x, &par[17]) + fitPeak(x, &par[20]) + fitPeak(x, &par[23]);
//	return fitBackground(x, &par[0]) \
//		+ fitPeak(x, &par[2]) + fitPeak(x, &par[5]) + fitPeak(x, &par[8]) + fitPeak(x, &par[11]) \
//		+ fitPeak2(x, &par[14], 0) + fitPeak2(x, &par[16], 1) + fitPeak2(x, &par[18], 2)+ fitPeak2(x, &par[20], 3);
//	return fitBackground(x, &par[0]) + fitPeak(x, &par[2]) + fitPeak(x, &par[5]) + fitPeak(x, &par[8]) + fitPeak(x, &par[11]);
}


// hist is a pointer to https://root.cern.ch/doc/master/classTH1D.html
void do_calibration(TH1F* hist, std::string pname, std::string pname_short, Double_t* peak_centres, Double_t* secondary_centres, Double_t* fitted_centres, Double_t &gain, Double_t &offset, Double_t &chisq, Double_t &ffchisq) {

	Double_t peak_constants[8];
	Double_t peak_means[8];
	Double_t peak_sigmas[8];

	fitted_centres[0]=0;
	fitted_centres[1]=0;
	fitted_centres[2]=0;
	fitted_centres[3]=0;
	
	gain=0;
	offset=0;
	chisq=0;
	ffchisq=0;
	
	std::string ver = "do_calibration_m2 v1";
	std::cout << ver << std::endl;

	std::cout << "Expected Max : " + std::to_string(hist->GetBinCenter(hist->GetMaximumBin())) << std::endl;
	std::cout << "The 4 peak centres : " << peak_centres[0] << "\t" << peak_centres[1] << "\t" << peak_centres[2] << "\t" << peak_centres[3] << std::endl;

//	// Set the legend
//	gStyle->SetOptFit(1);
//	gStyle->SetOptStat(1);

	// Set fit function options
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili");

	// Loop through the 4 main peak Q-value bin centres and fit Gaussian
	for (int i = 0; i < 4; i++) {

		// If we have a non-zero value we have a peak
		if (peak_centres[i] != 0) {

			TF1* peakfit = new TF1("peakfit", "gaus");
			hist->Fit("peakfit", "LQ+", "", peak_centres[i] - 2, peak_centres[i] + 2);

			peak_constants[i] = peakfit->GetParameter(0);  // y - height
			peak_means[i] = peakfit->GetParameter(1);      // x - centre
			peak_sigmas[i] = peakfit->GetParameter(2);     // width - sigma
		}
		else {
			peak_constants[i] = 0;
			peak_means[i] = 0;
			peak_sigmas[i] = 0;
		}
		std::cout << "FIT : " << peak_constants[i] << "\t" << peak_means[i] << "\t" << peak_sigmas[i] << std::endl;
	}

//	hist->Draw();

	//if (TRUE) { return; };

//	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili");

	std::cout << "Expected fine structure : " + std::to_string(hist->GetBinCenter(hist->GetMaximumBin())) << std::endl;
	std::cout << "The fine structure centres : " << secondary_centres[0] << "\t" << secondary_centres[1] << "\t" << secondary_centres[2] << "\t" << secondary_centres[3] << std::endl;

	double range = 0.9;
	// Loop through the fine structrure Q-value bin centres and fit Gaussian
	for (int i = 0; i < 4; i++) {

		// If we have a non-zero value we have a peak
		if (secondary_centres[i] != 0) {

			TF1* peakfit = new TF1("peakfit", "gaus");
			hist->Fit("peakfit", "LQ+", "", secondary_centres[i]-range-0.25, secondary_centres[i]+range);

			std::cout << "FIT* : Got " << peakfit->GetParameter(0) << ", " << peakfit->GetParameter(1) << ", " << peakfit->GetParameter(2);
			std::cout << " peak between " << secondary_centres[i] << " range " << range << std::endl;

			// Check how small the width is, if ok save the fit, if 'too small' (<2) just save the bin centre and set width to 1
			if (abs(secondary_centres[i] - peakfit->GetParameter(1)) < 2) {
				peak_constants[4 + i] = peakfit->GetParameter(0);  // y - height
				peak_means[4 + i] = peakfit->GetParameter(1);      // x - centre
				peak_sigmas[4 + i] = peakfit->GetParameter(2);     // width - sigma
			}
			else {
				peak_constants[4 + i] = hist->GetBinContent(secondary_centres[i]);  // y - height
				peak_means[4 + i] = secondary_centres[i];      // x - centre
				peak_sigmas[4 + i] = 1;     // width - sigma
			}

		}
		else {
			peak_constants[4 + i] = 0;
			peak_means[4 + i] = 0;
			peak_sigmas[4 + i] = 0;
		}
		//std::cout << "\nFIT* : " << peak_constants[4 + i] << "\t" << peak_means[4 + i] << "\t" << peak_sigmas[4 + i] << std::endl;
	}

	// Line below will draw graph of the fits of the Gaussians used by the search
//	hist->Draw();
//	return;

	// ********************************************
	// Graph the full fit
	// ********************************************

	// Set fit function options
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili");

	// create a TF1 with the range from 200 to 700 and n parameters
	TF1* fitFull = new TF1("fitFunction", fitFunction, 200, 700, 23);
	fitFull->SetNpx(500);
	fitFull->SetLineWidth(1);
	fitFull->SetLineColor(kMagenta);

	// Seed the background straight line fit
	fitFull->SetParameter(0, 1);
	fitFull->SetParameter(1, 1);

	// Set the Gaussian fit parameters for the 4 main peaks and first 3 fine structure peaks
	for (int i = 0; i < 7; i++) {
		std::cout << "** Adding params : " << 2 + 3 * i << " to " << 4 + 3 * i << std::endl;
		fitFull->SetParameter(2 + 3 * i, peak_constants[i]); // y - height
		fitFull->SetParameter(3 + 3 * i, peak_means[i]);     // x - centre
		fitFull->SetParameter(4 + 3 * i, peak_sigmas[i]);    // width - sigma
	}
// Separate processing of fine structure - no longer used as all 7 read in above loop
//	for (int i = 4; i < 4; i++) {
//		std::cout << "\n** Adding params : " << 14 + 2 * i << " to " << 15 + 2 * i << std::endl;
//		fitFull->SetParameter(14 + 2 * i, peak_constants[i]); // y - height
//		fitFull->SetParameter(15 + 2 * i, peak_means[i]);     // x - centre
//	}

	// Perform Full-fit
	hist->Fit("fitFunction", "Q", "E", 200, 700);
	ffchisq = fitFull->GetChisquare();
	cout << "Fullfit Reduced Chisq " << ffchisq << " / 23 " << std::endl;
	cout << "Fullfit Reduced Chisq " << ffchisq/23 << std::endl;

	// Guassian centres start at element 3
	fitted_centres[0] = fitFull->GetParameter(3+3*0);  // 1st centre
	fitted_centres[1] = fitFull->GetParameter(3+3*1);  // 2nd centre
	fitted_centres[2] = fitFull->GetParameter(3+3*2);  // 3rd centre
	fitted_centres[3] = fitFull->GetParameter(3+3*3);  // 4th centre	

	// Guassian centres start at element 4
	Double_t fitted_errors[4] = {0,0,0,0};
	fitted_errors[0] = fitFull->GetParameter(4+3*0);  // 1st error
	fitted_errors[1] = fitFull->GetParameter(4+3*1);  // 2nd error
	fitted_errors[2] = fitFull->GetParameter(4+3*2);  // 3rd error
	fitted_errors[3] = fitFull->GetParameter(4+3*3);  // 4th error	

	std::cout << "Full fit centres : " << fitted_centres[0] << "\t" << fitted_centres[1] << "\t" << fitted_centres[2] << "\t" << fitted_centres[3] << std::endl;
	std::cout << "Full fit errors : " << fitted_errors[0] << "\t" << fitted_errors[1] << "\t" << fitted_errors[2] << "\t" << fitted_errors[3] << std::endl;
	
	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=148Gd&unc=NDS
	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=239Pu&unc=NDS
	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=241Am&unc=NDS
	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=244Cm&unc=NDS
	// Peak Alpha Decay Q-Values (keV) for Gd-148, Pu-239, Am-241, Cm-244
	Double_t Q_alpha1[4] = { 3182.690, 5156.59, 5485.56, 5804.77 };
	Double_t Q_alpha1_errors[4] = { 0.01, 0.01, 0.01, 0.01 };

	// Set what's displayed in the legend
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(1);	
	hist->GetYaxis()->SetTitle("Frequency / counts");
	// Draw the full fit
	hist->Draw();

	// ********************************************
	// Graph of straigh line fit of Q-values vs ADC values
	// ********************************************

	// Set fit function options
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	
//	TGraph* gr1 = new TGraph(4, peak_centres, Q_alpha1);
	TGraph* gr1 = new TGraph(4, fitted_centres, Q_alpha1);
// Replace above with below if you want errors!
//	TGraph* gr1 = new TGraphErrors(4, fitted_centres, Q_alpha1, fitted_errors, Q_alpha1_errors);

	// !!!
	gr1->SetMarkerSize(1);
	gr1->SetMarkerStyle(88);

	// use "q" as 2nd argument to run quiet
	gr1->Fit("pol1", "Q");
	TF1* gr1_fit = (TF1*)gr1->GetFunction("pol1");
	// std::cout << gr1_fit->GetErrors() << std::endl;
	gain = gr1_fit->GetParameter(1);
	offset = gr1_fit->GetParameter(0); 	
	chisq = gr1_fit->GetChisquare();
	cout << "Offset/Gain Reduced Chisq " << chisq << " / 2 " << std::endl;
	cout << "Offset/Gain Reduced Chisq " << chisq/2 << std::endl;

	std::string title = pname + " - Peak Centroids";
	gr1->SetName(pname_short.data());
	gr1->SetTitle(title.data());
	gr1->GetXaxis()->SetTitle("ADC value");
	gr1->GetYaxis()->SetTitle("Q-Value / keV");
	gr1->GetXaxis()->SetLimits(400, 700);

	// Draw the gain & offset graph
	// Include next 2 lines calibration straight lines to draw fit graph
//	gr1->Draw();
//	gStyle->SetOptFit(1); // Display fit parameters in the legend

	std::cout << pname << " - Offset: " << offset << " Gain: " << gain << std::endl;

	return;
}

/*
 * Call with no parmeters to process all data for method 4 and write results file 
 * or specify 
 *	module_id, asic_id, and channel_id
 *	fit_method [1-4]
 *	write_results [0|1]
 */
void run_calibration(int module_id = -1, int asic_id = -1, int channel_id = -1, int fit_method = 0, int write_results = 1) {

	ofstream outfile;

	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=148Gd&unc=NDS
	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=239Pu&unc=NDS
	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=241Am&unc=NDS
	// https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc=244Cm&unc=NDS
	// Peak Alpha Decay Q-Values (keV) for Gd-148 (100%), Pu-239 (70.77%), Am-241 (84.80%), Cm-244 (76.90%)
	Double_t Q_alpha1[4] = { 3182.69, 5156.59, 5485.56, 5804.77 };

	// Other Alpha Decay Q-Values (keV) (over 10% intensity) 
	//  Gd-148 - None
	//  Pu-239 - 5144.3 (17.11%), 5105.5 (11.94%)
	//  Am-241 - 5442.80 (13.10%)
	//  Cm-244 - 5762.64 (23.10%) 
	Double_t Q_alpha2[4] = { 5144.3, 5442.80, 5762.64, 5105.5 };

	// used to record the full-fit centres
	Double_t fitted_centres[4] = { 0, 0, 0, 0 };
	Double_t gain;
	Double_t offset;
	Double_t chisq;
	Double_t ffchisq;

	int i_start; int i_end; int j_start; int j_end; int k_start; int k_end;

	// Chooose peak search method - set decsription, function pointer
	std::string fit_method_description;
	void (*search_func_ptr)(TH1F *, Double_t *, Double_t *, Double_t *, Double_t *) { nullptr };
	
	if (fit_method == 1) { // Read-ahead
		fit_method_description = "Method 1";
		search_func_ptr = &get_peak_centres_method1;
	}
	else if (fit_method == 2) { // Read-ahead with Gaussian
		fit_method_description = "Method 2";
		search_func_ptr = &get_peak_centres_method2;
	}
	else if (fit_method == 3) { // Start search based on Q Values
		fit_method_description = "Method 3";
		search_func_ptr = &get_peak_centres_method3;
	}
	else { // Start search based first Q-value and then known peak differences
		fit_method = 4;
		fit_method_description = "Method 4";
		search_func_ptr = &get_peak_centres_method4;
	}

	// Open file
//	TFile* f = new TFile("../data_files/R204_0.root", "READ");
	TFile* f = new TFile("../data_full/total_alpha_data.root", "READ");

	std::cout << "************* S T A R T *************" << std::endl;
	std::cout << "Using fit method " << fit_method_description << "(#" << fit_method << ")" << std::endl;

	// Declare pointers
	TH2F* hist2d;
	TH1F* hist1d;

	// Make Canvas
	TCanvas* c = new TCanvas();


	// If we have the default values of -1 run for all data  
	if (module_id == -1 && asic_id == -1 && channel_id == -1) {
		i_start = 0; i_end = 4; j_start = 0; j_end = 6; k_start = 0; k_end = 128;
		// Else run for specified module_id, asic_id, and channel_id  
	}
	else {
		i_start = module_id; i_end = module_id + 1; j_start = asic_id; j_end = asic_id + 1; k_start = channel_id; k_end = channel_id + 1;
	}

	if (write_results == 1) {
		outfile.open ("output.csv");
	}
	
	// Loop over modules
	for (int i = i_start; i < i_end; i++) {

		// Loop over asics
		for (int j = j_start; j < j_end; j++) {

			// Define histogram name
			std::string hname = "asic_hists/module_";
			hname += std::to_string(i);
			hname += "/asic_" + std::to_string(i);
			hname += "_" + std::to_string(j);

			// Get 2D histogram
			hist2d = f->Get<TH2F>(hname.data());

			// Loop over channel number
			for (int k = k_start; k < k_end; k++) {

				// Projection
//				std::string pname = "Projection_" + std::to_string(i);
//				pname += "_" + std::to_string(j);
//				pname += "_" + std::to_string(k);
//				pname += " - " + fit_method_description;

				// Projection
				std::string pname_short = "Channel " + std::to_string(k);
				pname_short += " - " + fit_method_description;

				std::string pname = "Module " + std::to_string(i);
				pname += " ASIC " + std::to_string(j);
				pname += " " + pname_short;

				// std::cout << hname + " " + pname << std::endl;

				hist1d = (TH1F*)hist2d->ProjectionY(pname_short.data(), k + 1, k + 1);

				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				std::cout << "Total Counts " << hist1d->GetEntries();

				hist1d->Draw();
				//hist1d->SetLineColor(kRed);`


				// Do calibration
				Double_t peak_centres[4];
				Double_t secondary_centres[4];

				// Search for peaks using the search function of the chosen method
				search_func_ptr(hist1d, Q_alpha1, Q_alpha2, peak_centres, secondary_centres);
				
				// Perform calibration using peaks, return fitted centres, gan and offset
				do_calibration(hist1d, pname.data(), pname_short.data(), peak_centres, secondary_centres, fitted_centres, gain, offset, chisq, ffchisq);

				if (write_results == 1) {
					outfile << i << "," << j << "," << k << "," << hist1d->GetEntries() << "," << gain << "," << offset << "," << chisq << "," << ffchisq << ",";
					outfile << fitted_centres[0] << "," << fitted_centres[1] << "," << fitted_centres[2] << "," << fitted_centres[3];
					outfile <<  std::endl;
				}
			} // k

		} // j

	} // i


	if (write_results == 1) {
		outfile.close();
	}

	// Close file
	//f->Close();

	return;
}


//
// Useful Links
//
// https://root.cern/doc/master/pyroot_2multifit_8py.html
