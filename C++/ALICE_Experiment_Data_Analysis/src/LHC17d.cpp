double CB22(double *x, double *par)
{
	double m = par[1];
	double s = par[2];
	double n = par[3];
	double a = par[4];
	double kMean = 0.135;
	double dx = (x[0] - m) / s;
	if (dx > -a)
	{
		return par[0] * exp(-dx * dx / 2.) +
		par[5] +
		par[6] * (x[0] - kMean) +
		par[7] * (x[0] - kMean) * (x[0] - kMean);
	}
	else
	{
		double A = TMath::Power((n / TMath::Abs(a)), n) * TMath::Exp(-a * a / 2);
		double B = n / TMath::Abs(a) - TMath::Abs(a);
		return par[0] * A * TMath::Power((B - dx), -n) +
		par[5] +
		par[6] * (x[0] - kMean) +
		par[7] * (x[0] - kMean) * (x[0] - kMean);
	}
}

double pol22(double *x, double *par)
{

	double kMean = 0.135;
	return par[0] +
		   par[1] * (x[0] - kMean) +
		   par[2] * (x[0] - kMean) * (x[0] - kMean);
}

void LHC17d(){

	gStyle->SetOptStat(0);

	TFile *file1 = new TFile("LHC17d.root", "read");

	THashList *L1 = (THashList*)file1->Get("TaggingkINT7V0M");

	TH2F *hist1 = (TH2F*)L1->FindObject("hInvM_Re_Emin1_All_cent0");


	TH2F *CloneReal = (TH2F*)hist1->Clone();
	TH2F *CloneResult = (TH2F*)hist1->Clone();

	TH2F *CloneRealGaus = (TH2F*)hist1->Clone();
	TH2F *CloneResultGaus = (TH2F*)hist1->Clone();

	TH2F *hist2 = (TH2F*)L1->FindObject("hInvM_Mi_Emin1_All_cent0");


	TCanvas *canvas1 = new TCanvas("canvas1", " ");
	canvas1->cd();
	hist1->SetTitle("Two-photon inv. mass vs pion P_{T}");
	hist1->GetYaxis()->SetRangeUser(0, 20);
	hist1->GetYaxis()->SetTitle("P_{T}, GeV");
	hist1->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV");
	hist1->Draw("colz");

	canvas1->cd();

	TPaveText *pave = new TPaveText(0.2, 0.7, 0.6, 0.8, "NDC");
	pave->SetFillColor(0); 
	pave->SetLineColor(kRed); 
	pave->SetLineWidth(2); 
	pave->AddText("p+p->#pi^{0} + X, #sqrt{S} = 13 TeV");
	pave->SetTextColor(kRed);
	pave->Draw();


	TCanvas *canvas2 = new TCanvas("canvas2", " ");
	canvas2->cd();
	hist2->SetTitle("Two-photon inv. mass vs pion P_{T}");
	hist2->GetYaxis()->SetRangeUser(0, 20);
	hist2->GetYaxis()->SetTitle("P_{T}, GeV");
	hist2->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV");
	hist2->Draw("colz");

	canvas2->cd();
	pave->Draw();

	TH2F *hist3 = (TH2F*)L1->FindObject("hSelEvents");

	TH2F *CloneMixed = (TH2F*)hist2->Clone();
	TH2F *CloneNormedMixed = (TH2F*)hist2->Clone();

	hist1->Divide(hist2); 

	TH1D *Pt = hist1->ProjectionY();
	Int_t nbins = Pt->GetNbinsX();

	// double low_pt [] = {1.0, 1.5, 2.0, 3.0, 4.0, 4.5};
	// double up_pt [] =  {1.5, 2.0, 3.0, 4.0, 4.5, 10.0};

	double low_pt [] = {1.0, 1.4, 2.0, 2.5, 3.2, 5.0};
	double up_pt [] =  {1.4, 2.0, 2.5, 3.2, 5.0, 10};	

	// double low_pt [] = {0.3, 0.8, 1.3, 1.8, 2.3, 3.1};
	// double up_pt [] =  {0.8, 1.3, 1.8, 2.3, 3.1, 4.0};


	int number_of_proj = sizeof(low_pt)/ sizeof(low_pt[0]); cout << number_of_proj << endl;

	double low_edge[nbins];
	double up_edge[nbins];

	for (int i = 1; i <= nbins; ++i)
	{
		low_edge[i] = hist1->ProjectionY()->GetBinLowEdge(i);
		up_edge[i] = hist1->ProjectionY()->GetBinLowEdge(i) + hist1->ProjectionY()->GetBinWidth(i);
	}

	double mid_pt [number_of_proj];
	double mid_pt_er[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		mid_pt[i] = (low_pt[i] + up_pt[i])/2;
		mid_pt_er[i] = 0;
	}

	double low_pt_bin[number_of_proj];
	double up_pt_bin[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		low_pt_bin[i] = hist1->GetYaxis()->FindBin(low_pt[i]);
		up_pt_bin[i] = hist1->GetYaxis()->FindBin(up_pt[i]-0.1*up_edge[i]);
		cout << low_pt_bin[i] << "\t" << up_pt_bin[i] << endl;
	}

	TH1D *mixed[number_of_proj];
	TH1D *real[number_of_proj];
	TH1D *real_mixed[number_of_proj];

	double Norm1[number_of_proj];
	double Norm2[number_of_proj];
	double Norm3[number_of_proj];

	double number[number_of_proj];
	double peaks[number_of_proj];
	double sigma[number_of_proj];

	double peaks_er[number_of_proj];
	double sigma_er[number_of_proj];

	int j = 3;	
	TCanvas *c[j];
	for (int i = 0; i < j; ++i)
	{
		c[i] = new TCanvas(Form("c%d", i), Form("c%d", i));
		c[i]->Divide(3, 2);
	}

	double max_bin[number_of_proj];
	double chisqure[number_of_proj];
	double NDF[number_of_proj];

	TF1 *func[number_of_proj];
	TF1 *BG[number_of_proj];

	double board1 = 0.06;
	double board2 = 0.40;

	// double board1 = 0.08;
	// double board2 = 0.35;

	for (int i = 0; i < number_of_proj; ++i)
	{
		real[i] = CloneReal->ProjectionX(Form("real_projections(No Field)%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		real[i] -> SetTitle(Form("real to normed_mixed (Crystalball, No Field) %.1f < Pt < %.1f", low_pt[i], up_pt[i]));
		// c[0]->cd(i+1);
		// c[0]->SetGrid();
		// real[i]->Draw("same");

		mixed[i] = CloneMixed->ProjectionX(Form("mixed_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
	 	mixed[i] -> SetTitle(Form("mixed distribution %.1f < Pt < %.1f", low_pt[i], up_pt[i]));
	 	// c[1] -> cd(i+1);
	 	// mixed[i]->Draw("same");

	 	// c[2]->cd(i+1);
		// mixed[i]->Draw("same");
		// real[i]->Draw("same");

	 	real_mixed[i] = CloneReal->ProjectionX(Form("real/mixed_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		real_mixed[i] -> SetTitle(Form("%.1f < Pt < %.1f", low_pt[i], up_pt[i])); 
		real_mixed[i]->Divide(mixed[i]);
		real_mixed[i]->GetXaxis()->SetRangeUser(0.06, 0.40);
		c[0]->cd(i+1);
		real_mixed[i]->Draw("same");

		func[i] = new TF1("func", CB22, 0., 1.0, 8);
		BG[i] = new TF1("BG", "[0] + [1]*x + [2]*x*x", 0.04, 1);


		for (int k = 0; k < real_mixed[i]->GetNbinsX(); k++)
		{
			if ( k > 50 && k < 90)
			{
				
				Double_t test=real_mixed[i]->GetBinContent(k);
				if(test > max_bin[i]){
					max_bin[i] = test;
			   }
			}
		}


		if (i == 0)
		{
			func[i]->SetParameters(max_bin[i]*2,
			0.135,
			0.01,
			0.05,
			1.7,
			max_bin[i] / 40.,
			-max_bin[i]/ 10.,
			max_bin[i]/ 10.);
		}
		else if (i == 1)
		{
			func[i]->SetParameters(max_bin[i]*2,
			0.135,
			0.006,
			0.05,
			1.7,
			max_bin[i] / 40.,
			-max_bin[i]/ 10.,
			max_bin[i]/ 10.);	
		}

		else if (i > 1 && i != 5)
		{
			func[i]->SetParameters(max_bin[i],
			0.135,
			0.006,
			0.05,
			1.7,
			max_bin[i] / 40.,
			-max_bin[i]/ 10.,
			max_bin[i]/ 10.);
		}

		else if (i == 5)
		{
			func[i]->SetParameters(max_bin[i]/0.8,
			0.135,
			0.006,
			0.05,
			1.7,
			max_bin[i] / 40.,
			-max_bin[i]/ 10.,
			max_bin[i]/ 10.);
		}


		// if (i == 0)
		// {
		// 	func[i]->SetParameters(max_bin[i]*2,
		// 	0.135,
		// 	0.01,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }
		// else if (i == 1)
		// {
		// 	func[i]->SetParameters(max_bin[i]*2,
		// 	0.135,
		// 	0.006,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);	
		// }

		// else if (i > 1)
		// {
		// 	func[i]->SetParameters(max_bin[i],
		// 	0.135,
		// 	0.006,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }
		// else if ( i == 5)
		// {
		// 	func[i]->SetParameters(max_bin[i]/2,
		// 	0.135,
		// 	0.006,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }

		// if (i == 0)
		// {
		// 	func[i]->SetParameters(max_bin[i]/1.1,
		// 	0.135,
		// 	0.01,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);

		// 	BG[i]->SetParameters(max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }
		// else if (i == 1)
		// {
		// 	func[i]->SetParameters(max_bin[i]/1.1,
		// 	0.135,
		// 	0.01,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);

		// 	BG[i]->SetParameters(max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);	
		// }
		// else if(i == 2){
		// 	func[i]->SetParameters(max_bin[i]*2,
		// 	0.135,
		// 	0.006,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);

		// 	BG[i]->SetParameters(max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }
		// else if (i > 2)
		// {
		// 	func[i]->SetParameters(max_bin[i],
		// 	0.135,
		// 	0.006,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);

		// 	BG[i]->SetParameters(max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }	

		cout << "**********************" << endl;
		cout << i+1  << " -я проекция " << low_pt[i] << " Pt < " << up_pt[i] << endl;
		cout << "**********************" << endl;

		// if (i != 5)
		// {
		// 	real_mixed[i]->Fit("func", "", "", board1, board2);
		// }
		// else if (i == 5)
		// {
		// 	real_mixed[i]->Fit("func", "", "", 0.1, 0.25);
		// }
		real_mixed[i]->Fit("func", "", "", board1, board2); // было 0.06
		BG[i]->SetLineColor(kViolet);
		//real_mixed[i]->Fit("BG", "R+");

		chisqure[i] = func[i]->GetChisquare();
		NDF[i] = func[i]->GetNDF();

		Norm1[i] = func[i]->GetParameter(5);
		Norm2[i] = func[i]->GetParameter(6);
		Norm3[i] = func[i]->GetParameter(7);  
	}	

	for (int i = 0; i < number_of_proj; ++i)
	{
		cout << "Chisquare/NDF = " << chisqure[i] << "/" << NDF[i] << " = " << chisqure[i]/NDF[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	}

	TF1 *NormFunc[number_of_proj];
	TF1 *ResultFunc[number_of_proj];
	TH1D *NormMixed[number_of_proj];
	TH1D *result[number_of_proj];

	double max_bin1[number_of_proj];
	double chisqure1[number_of_proj];
	double NDF1[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		NormMixed[i] = CloneNormedMixed->ProjectionX(Form("NormMixed_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		NormMixed[i] -> SetTitle(Form("normed_mixed distribution %.1f < Pt < %.1f", low_pt[i], up_pt[i]));

		NormFunc[i] = new TF1("NormFunc", "[0] + [1]*(x - 0.135) +[2]*(x-0.135)*(x-0.135)" , 0., 1.);

		NormFunc[i]->SetParameter(0, Norm1[i]);
		NormFunc[i]->SetParameter(1, Norm2[i]);
		NormFunc[i]->SetParameter(2, Norm3[i]);

		NormMixed[i]->Multiply(NormFunc[i]);

		// c[1]->cd(i+1);
		// NormMixed[i]->Draw("same");


		result[i] = CloneResult->ProjectionX(Form("result_projections(No Field)%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		result[i] -> SetTitle(Form(" %.1f < Pt < %.1f", low_pt[i], up_pt[i]));

		result[i]->Add(NormMixed[i], -1);

		c[1]->cd(i+1);
		c[1]->SetGrid(i+1);

		real[i]->SetAxisRange(board1, board2, "X");
		NormMixed[i]->SetAxisRange(board1, board2,"X");

		real[i]->Draw("same");
		NormMixed[i]->Draw("same");

		ResultFunc[i] =  new TF1("ResultFunc", CB22, 0., 1.0, 8);

		for (int k = 0; k < result[i]->GetNbinsX(); k++)
		{
			if ( k > 50 && k < 90)
			{
				
				Double_t test=result[i]->GetBinContent(k);
				if(test > max_bin1[i]){
					max_bin1[i] = test;
			   }
			}
		}

		c[2]->cd(i+1);

		ResultFunc[i]->SetParameters(max_bin1[i],
		0.135,
		0.006,
		0.05,
		1.7,
		0,
		0,
		0);

		result[i]->Fit("ResultFunc", "", "", board1, board2);

		c[2]->SetGrid();
		result[i]->SetAxisRange(board1, board2,"X");

		result[i]->Draw("same");

		chisqure1[i] = ResultFunc[i]->GetChisquare();
		NDF1[i] = ResultFunc[i]->GetNDF();

		peaks[i] = ResultFunc[i]->GetParameter(1);
		sigma[i] = ResultFunc[i]->GetParameter(2);

		peaks_er[i] = ResultFunc[i]->GetParError(1);
   		sigma_er[i] = ResultFunc[i]->GetParError(2);
 
	}


	for (int i = 0; i < number_of_proj; ++i)
	{
		cout << "Chisquare1/NDF1 = " << chisqure1[i] << "/" << NDF1[i] << " = " << chisqure1[i]/NDF1[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	}

	// TFile *ResultFile = new TFile("compare_low_Pt.root", "update");
	// for (int i = 0; i < number_of_proj; ++i)
	// {
	// 	real[i]->Divide(mixed[i]);
	// 	real[i]->Write();
	// 	result[i]->Write();
	// }
    // ResultFile->Close();



    //*****Method1. Hist integral****
	// double number_of_pairs[number_of_proj];
	// double number_of_pairs_er[number_of_proj];

	// for (int i = 0; i < number_of_proj; ++i)
	// {
	// 	for (int k = result[i]->FindBin(peaks[i] - 2.5*sigma[i]); k < result[i]->FindBin(peaks[i] + 2.5*sigma[i]); ++k)
	// 	{
	// 		number_of_pairs[i] += result[i]->GetBinContent(k)*0.002;
	// 		number_of_pairs_er[i] += result[i]->GetBinError(k)*0.002;
	// 	}
	// }
	// //*****Method1. Hist integral****


	//*****Method1. Hist integral****
	double number_of_pairs[number_of_proj];
	double epsilon[number_of_proj];
	double number_of_pairs_er[number_of_proj];
	double number_of_pairs_er22[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		if (i != 1 && i != 2 && i!= 3 && i!= 4 && i != 5)
		{
			for (int k = result[i]->FindBin(peaks[i] - 2.5*sigma[i]); k < result[i]->FindBin(peaks[i] + 2.5*sigma[i]); ++k)
			{
				number_of_pairs[i] += result[i]->GetBinContent(k)*0.002;
				epsilon[i] += 1/sqrt(result[i]->GetBinContent(k));
				number_of_pairs_er[i] = epsilon[i]*number_of_pairs[i]*0.002;
				number_of_pairs_er22[i] += result[i]->GetBinError(k)*0.002;
			}
		}
		else if(i == 1)
		{
			for (int k = result[i]->FindBin(peaks[i] - 1.1*sigma[i]); k < result[i]->FindBin(peaks[i] + 1.1*sigma[i]); ++k)
			{
				number_of_pairs[i] += result[i]->GetBinContent(k)*0.002;
				epsilon[i] += 1/sqrt(result[i]->GetBinContent(k));
				number_of_pairs_er[i] = epsilon[i]*number_of_pairs[i]*0.002;
				number_of_pairs_er22[i] += result[i]->GetBinError(k)*0.002;
			}
		}
		else if(i == 2)
		{
			for (int k = result[i]->FindBin(peaks[i] - 4.0*sigma[i]); k < result[i]->FindBin(peaks[i] + 4.0*sigma[i]); ++k)
			{
				number_of_pairs[i] += result[i]->GetBinContent(k)*0.002;
				epsilon[i] += 1/sqrt(result[i]->GetBinContent(k));
				number_of_pairs_er[i] = epsilon[i]*number_of_pairs[i]*0.002;
				number_of_pairs_er22[i] += result[i]->GetBinError(k)*0.002;
			}
		}
		else if(i == 3)
		{
			for (int k = result[i]->FindBin(peaks[i] - 1.6*sigma[i]); k < result[i]->FindBin(peaks[i] + 1.6*sigma[i]); ++k)
			{
				number_of_pairs[i] += result[i]->GetBinContent(k)*0.002;
				epsilon[i] += 1/sqrt(result[i]->GetBinContent(k));
				number_of_pairs_er[i] = epsilon[i]*number_of_pairs[i]*0.002;
				number_of_pairs_er22[i] += result[i]->GetBinError(k)*0.002;
			}
		}
		else if(i == 4)
		{
			for (int k = result[i]->FindBin(peaks[i] - 1.9*sigma[i]); k < result[i]->FindBin(peaks[i] + 1.9*sigma[i]); ++k)
			{
				number_of_pairs[i] += result[i]->GetBinContent(k)*0.002;
				epsilon[i] += 1/sqrt(result[i]->GetBinContent(k));
				number_of_pairs_er[i] = epsilon[i]*number_of_pairs[i]*0.002;
				number_of_pairs_er22[i] += result[i]->GetBinError(k)*0.002;
			}
		}
		else if (i == 5)
		{
			for (int k = result[i]->FindBin(peaks[i] - 3.9*sigma[i]); k < result[i]->FindBin(peaks[i] + 3.9*sigma[i]); ++k)
			{
				number_of_pairs[i] += result[i]->GetBinContent(k)*0.002;
				epsilon[i] += 1/sqrt(result[i]->GetBinContent(k));
				number_of_pairs_er[i] = epsilon[i]*number_of_pairs[i]*0.002;
				number_of_pairs_er22[i] += result[i]->GetBinError(k)*0.002;
			}
		}
	}
	//*****Method1. Hist integral****

	cout << "&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	for (int i = 0; i < number_of_proj; ++i)
	{
		cout << epsilon[i] << "\t" << number_of_pairs_er[i] << endl;
	}
	cout << "&&&&&&&&&&&&&&&&&&&&&&&" << endl;


	//*****Method2. Function integral****
	double number_of_pairs2[number_of_proj];
	double number_of_pairs2_er[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		number_of_pairs2[i] = ResultFunc[i]->Integral(peaks[i] - 3*sigma[i], peaks[i] + 3*sigma[i]);
		//number_of_pairs2_er[i] = ResultFunc[i]->IntegralError(peaks[i] - 3*sigma[i], peaks[i] + 3*sigma[i]);
		//number_of_pairs2_er[i] = number_of_pairs_er[i] - 0.03*number_of_pairs_er[i];
		number_of_pairs2_er[i] = 0.05*number_of_pairs2[i];

	}
	//*****Method2. Function integral****











	//**************теперь то же самое гауссом***********************


	TH2F *CloneMixedGaus = (TH2F*)hist2->Clone();
	TH2F *CloneNormedMixedGaus = (TH2F*)hist2->Clone();

	TH1D *mixedGaus[number_of_proj];
	TH1D *realGaus[number_of_proj];
	TH1D *real_mixedGaus[number_of_proj];

	double Norm1Gaus[number_of_proj];
	double Norm2Gaus[number_of_proj];
	double Norm3Gaus[number_of_proj];

	double peaksGaus[number_of_proj];
	double sigmaGaus[number_of_proj];

	double peaks_erGaus[number_of_proj];
	double sigma_erGaus[number_of_proj];

	int jGaus = 3;	
	TCanvas *cGaus[jGaus];
	for (int i = 0; i < jGaus; ++i)
	{
		cGaus[i] = new TCanvas(Form("cGaus%d", i), Form("cGaus%d", i));
		cGaus[i]->Divide(3, 2);
	}

	double max_binGaus[number_of_proj];
	double chisqureGaus[number_of_proj];
	double NDFGaus[number_of_proj];

	TF1 *funcGaus[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		realGaus[i] = CloneRealGaus->ProjectionX(Form("realGaus_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		realGaus[i] -> SetTitle(Form("real to normed_mixed (Gaus) %.1f < Pt < %.1f", low_pt[i], up_pt[i]));

		realGaus[i]->SetAxisRange(0.06, 0.35,"X");

		mixedGaus[i] = CloneMixedGaus->ProjectionX(Form("mixedGaus_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
	 	mixedGaus[i] -> SetTitle(Form("mixedGaus distribution %.1f < Pt < %.1f", low_pt[i], up_pt[i]));

	 	real_mixedGaus[i] = CloneRealGaus->ProjectionX(Form("real/mixedGaus_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		real_mixedGaus[i] -> SetTitle(Form("real/mixedGaus distribution %.1f < Pt < %.1f", low_pt[i], up_pt[i])); 
		real_mixedGaus[i]->Divide(mixedGaus[i]);
		cGaus[0]->cd(i+1);
		real_mixedGaus[i]->Draw("same");


		funcGaus[i] = new TF1("funcGaus", "gaus(0) + [3] + [4]*(x - 0.135) +[5]*(x-0.135)*(x-0.135)", 0., 1.0);

		for (int k = 0; k < real_mixedGaus[i]->GetNbinsX(); k++)
		{
			if ( k > 50 && k < 90)
			{
				
				Double_t test=real_mixedGaus[i]->GetBinContent(k);
				if(test > max_binGaus[i]){
					max_binGaus[i] = test;
			   }
			}
		}

		funcGaus[i]->SetParameter(0, max_binGaus[i]);
		funcGaus[i]->SetParameter(1, 0.135);
		funcGaus[i]->SetParameter(2, 0.006);
		funcGaus[i]->SetParameter(3, max_bin[i] / 40.);
		funcGaus[i]->SetParameter(4, -max_bin[i]/ 10.);
		funcGaus[i]->SetParameter(5, max_bin[i]/ 10.);

		cout << "**********************" << endl;
		cout << i+1  << " -я проекция " << low_pt[i] << " Pt < " << up_pt[i] << endl;
		cout << "**********************" << endl;

		real_mixedGaus[i]->Fit("funcGaus", "", "", 0.06, 0.35);

		chisqureGaus[i] = funcGaus[i]->GetChisquare();
		NDFGaus[i] = funcGaus[i]->GetNDF();

		Norm1Gaus[i] = funcGaus[i]->GetParameter(3);
		Norm2Gaus[i] = funcGaus[i]->GetParameter(4);
		Norm3Gaus[i] = funcGaus[i]->GetParameter(5);  
	}	

	for (int i = 0; i < number_of_proj; ++i)
	{
		cout << "ChisquareGaus/NDFGaus = " << chisqureGaus[i] << "/" << NDFGaus[i] << " = " << chisqureGaus[i]/NDFGaus[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	}

	TF1 *NormFuncGaus[number_of_proj];
	TF1 *ResultFuncGaus[number_of_proj];
	TH1D *NormMixedGaus[number_of_proj];
	TH1D *resultGaus[number_of_proj];

	double max_bin1Gaus[number_of_proj];
	double chisqure1Gaus[number_of_proj];
	double NDF1Gaus[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		NormMixedGaus[i] = CloneNormedMixedGaus->ProjectionX(Form("NormMixedGaus_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		NormMixedGaus[i] -> SetTitle(Form("normed_mixedGaus distribution %.1f < Pt < %.1f", low_pt[i], up_pt[i]));

		NormFuncGaus[i] = new TF1("NormFuncGaus", "[0] + [1]*(x - 0.135) +[2]*(x-0.135)*(x-0.135)" , 0., 1.);

		NormFuncGaus[i]->SetParameter(0, Norm1Gaus[i]);
		NormFuncGaus[i]->SetParameter(1, Norm2Gaus[i]);
		NormFuncGaus[i]->SetParameter(2, Norm3Gaus[i]);

		NormMixedGaus[i]->Multiply(NormFuncGaus[i]);

		resultGaus[i] = CloneResultGaus->ProjectionX(Form("resultGaus_projections%d", i+1), low_pt_bin[i], up_pt_bin[i]);
		resultGaus[i] -> SetTitle(Form("resultGaus distribution %.1f < Pt < %.1f", low_pt[i], up_pt[i]));

		resultGaus[i]->Add(NormMixedGaus[i], -1);

		cGaus[1]->cd(i+1);


		realGaus[i]->Draw("same");
		NormMixedGaus[i]->Draw("same");

		ResultFuncGaus[i] =  new TF1("ResultFuncGaus", "gaus(0)", 0., 1.0);

		for (int k = 0; k < resultGaus[i]->GetNbinsX(); k++)
		{
			if ( k > 50 && k < 90)
			{
				
				Double_t test=result[i]->GetBinContent(k);
				if(test > max_bin1Gaus[i]){
					max_bin1Gaus[i] = test;
			   }
			}
		}

		ResultFuncGaus[i]->SetParameter(0, max_bin1Gaus[i]);
		ResultFuncGaus[i]->SetParameter(1, 0.135);
		ResultFuncGaus[i]->SetParameter(2, 0.006);

		cGaus[2]->cd(i+1);

		resultGaus[i]->Fit("ResultFuncGaus", "", "", 0.06, 0.35);

		cGaus[2]->SetGrid();
		resultGaus[i]->SetAxisRange(0.06, 0.35,"X");
		NormMixedGaus[i]->SetAxisRange(0.06, 0.35,"X");

		resultGaus[i]->Draw("same");

		chisqure1Gaus[i] = ResultFuncGaus[i]->GetChisquare();
		NDF1Gaus[i] = ResultFuncGaus[i]->GetNDF();

		peaksGaus[i] = ResultFuncGaus[i]->GetParameter(1);
		sigmaGaus[i] = ResultFuncGaus[i]->GetParameter(2);

		peaks_erGaus[i] = ResultFuncGaus[i]->GetParError(1);
   		sigma_erGaus[i] = ResultFuncGaus[i]->GetParError(2);
	}

	for (int i = 0; i < number_of_proj; ++i)
	{
		cout << "Chisquare1Gaus/NDF1Gaus = " << chisqure1Gaus[i] << "/" << NDF1Gaus[i] << " = " << chisqure1Gaus[i]/NDF1Gaus[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	}

	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7" << endl;

	// cout << "***********************************" << endl;
	// for (int i = 0; i < number_of_proj; ++i)
	// {
	// 	cout << "Chisquare/NDF = " << chisqure[i] << "/" << NDF[i] << " = " << chisqure[i]/NDF[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	// }
	// cout << "***********************************" << endl;

	// cout << "***********************************" << endl;
	// for (int i = 0; i < number_of_proj; ++i)
	// {
	// 	cout << "Chisquare1/NDF1 = " << chisqure1[i] << "/" << NDF1[i] << " = " << chisqure1[i]/NDF1[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	// }
	// cout << "***********************************" << endl;

	// cout << "***********************************" << endl;
	// for (int i = 0; i < number_of_proj; ++i)
	// {
	// 	cout << "ChisquareGaus/NDFGaus = " << chisqureGaus[i] << "/" << NDFGaus[i] << " = " << chisqureGaus[i]/NDFGaus[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	// }
	// cout << "***********************************" << endl;

	// cout << "***********************************" << endl;
	// for (int i = 0; i < number_of_proj; ++i)
	// {
	// 	cout << "Chisquare1Gaus/NDF1Gaus = " << chisqure1Gaus[i] << "/" << NDF1Gaus[i] << " = " << chisqure1Gaus[i]/NDF1Gaus[i] << "\t" << low_pt[i] << " < Pt < " << up_pt[i] << endl;
	// }
	// cout << "***********************************" << endl;

	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7" << endl;

	//*****Method1. Hist integral****
	double number_of_pairsGaus[number_of_proj];
	double number_of_pairs_erGaus[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		for (int k = resultGaus[i]->FindBin(peaksGaus[i] - 3*sigmaGaus[i]); k < resultGaus[i]->FindBin(peaksGaus[i] + 3*sigmaGaus[i]); ++k)
		{
			number_of_pairsGaus[i] += resultGaus[i]->GetBinContent(k)*0.002;
			number_of_pairs_erGaus[i] += resultGaus[i]->GetBinError(k)*0.002;
		}
	}
	//*****Method1. Hist integral****


	//*****Method2. Function integral****
	double number_of_pairs2Gaus[number_of_proj];
	double number_of_pairs2_erGaus[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		number_of_pairs2Gaus[i] = ResultFunc[i]->Integral(peaksGaus[i] - 3*sigmaGaus[i], peaksGaus[i] + 3*sigmaGaus[i]);
		//number_of_pairs2_er[i] = ResultFunc[i]->IntegralError(peaks[i] - 3*sigma[i], peaks[i] + 3*sigma[i]);
		//number_of_pairs2_erGaus[i] = number_of_pairs_erGaus[i] - 0.03*number_of_pairs_erGaus[i];
		number_of_pairs2_erGaus[i] = 0.05*number_of_pairs2Gaus[i];
	}
	//*****Method2. Function integral****

	// for (int i = 0; i < hist3->GetNbinsX(); ++i)
	// {
	// 	cout << hist3->GetBinContent(i+1) << "\t" << hist3->GetBinLowEdge(i+1) << endl;
	// }

	double normed_koeff = hist3->GetBinContent(7);

	double normed_number_of_pairs[number_of_proj];
	double normed_number_of_pairs2[number_of_proj];
	double normed_number_of_pairsGaus[number_of_proj];
	double normed_number_of_pairs2Gaus[number_of_proj];

	double normed_number_of_pairs_er[number_of_proj];
	double normed_number_of_pairs2_er[number_of_proj];
	double normed_number_of_pairs_erGaus[number_of_proj];
	double normed_number_of_pairs_er2Gaus[number_of_proj];

	for (int i = 0; i < number_of_proj; ++i)
	{
		normed_number_of_pairs[i] = number_of_pairs[i]/normed_koeff;
		normed_number_of_pairs2[i] = number_of_pairs2[i]/normed_koeff;
		normed_number_of_pairsGaus[i] = number_of_pairsGaus[i]/normed_koeff;
		normed_number_of_pairs2Gaus[i] = number_of_pairs2Gaus[i]/normed_koeff;

		normed_number_of_pairs_er[i] = number_of_pairs_er[i]/normed_koeff;
		normed_number_of_pairs2_er[i] = number_of_pairs2_er[i]/normed_koeff;
		normed_number_of_pairs_erGaus[i] = number_of_pairs_erGaus[i]/normed_koeff;
		normed_number_of_pairs_er2Gaus[i] = number_of_pairs2_erGaus[i]/normed_koeff;
	}

	int number_of_graphs = 8;

	TGraphErrors *gr[number_of_graphs];
	gr[0] = new TGraphErrors(number_of_proj, mid_pt, peaks, mid_pt_er, peaks_er);
	gr[1] = new TGraphErrors(number_of_proj, mid_pt, sigma, mid_pt_er,  sigma_er);
	gr[2] = new TGraphErrors(number_of_proj, mid_pt, normed_number_of_pairs, mid_pt_er,  normed_number_of_pairs_er);
	gr[3] = new TGraphErrors(number_of_proj, mid_pt, normed_number_of_pairs2, mid_pt_er,  normed_number_of_pairs2_er);

	gr[4] = new TGraphErrors(number_of_proj, mid_pt, peaksGaus, mid_pt_er, peaks_erGaus);
	gr[5] = new TGraphErrors(number_of_proj, mid_pt, sigmaGaus, mid_pt_er,  sigma_erGaus);
	gr[6] = new TGraphErrors(number_of_proj, mid_pt, normed_number_of_pairsGaus, mid_pt_er,  normed_number_of_pairs_erGaus);
	gr[7] = new TGraphErrors(number_of_proj, mid_pt, normed_number_of_pairs2Gaus, mid_pt_er,  normed_number_of_pairs_er2Gaus);

	string name [] = {"m_{#gamma#gamma}", "#sigma", "number of pairs N in peaks (hist integral)", "number of pairs N in peaks (func integral)", "m_{#gamma#gamma}", "#sigma", "number of pairs N in peaks (hist integral)", "number of pairs N in peaks (func integral)"};
	string name2 [] = {", GeV", ", GeV", " ", " ", ", GeV", ", GeV", " ", " "};
	string name3 [] = {"m_{#gamma#gamma}", "#sigma", " N ", " N ", "m_{#gamma#gamma}", "#sigma", " N ", " N "};

	for (int i = 0; i < number_of_graphs; ++i)
	{
		gr[i]->SetTitle(("Dependency graph " + name[i] + " from P_{T}").c_str());
		gr[i]->GetXaxis()->SetTitle("P_{T}, GeV");

		gr[i]->GetYaxis()->SetTitle((name3[i] + name2[i]).c_str());

		if (i > 3)
		{
			gr[i]->SetMarkerStyle(30);
			gr[i]->SetMarkerSize(1.5);
			gr[i]->SetMarkerColor(3);
			gr[i]->SetLineColor(3);
			gr[i]->SetLineStyle(0);
		}
		else{
			gr[i]->SetMarkerStyle(45);
			gr[i]->SetMarkerSize(1.5);
			gr[i]->SetMarkerColor(6);
			gr[i]->SetLineColor(6);
			gr[i]->SetLineStyle(0);
		}

	}

	// for (int i = 0; i < number_of_proj; ++i)
	// {
	// 	cout << sigmaGaus[i] << "\t" << sigma[i] << endl;
	// }

	//**********compare**********

	int num_of_compare = 4;
	TCanvas *compare[num_of_compare];

	for (int i = 0; i < num_of_compare; ++i)
	{
		compare[i] = new TCanvas(Form("compare%d", i + 1), Form("compare%d", i+1));
	}

	compare[0]->cd();
	compare[0]->SetGrid();
	gr[0]->Draw("AP");
	gr[4]->Draw("P");

	TLegend *l1 = new TLegend(0.6, 0.6, 0.75, 0.75);
    l1->SetHeader("compare1 ","C"); 
    l1->AddEntry(gr[0],"Crystalball", "lep");
    l1->AddEntry(gr[4],"gaus", "lp");
    l1->Draw();

	compare[1]->cd();
	compare[1]->SetGrid();
	gr[1]->Draw("AP");
	gr[5]->Draw("P");


	TLegend *l2 = new TLegend(0.6, 0.6, 0.75, 0.75);
    l2->SetHeader("compare2 ","C"); 
    l2->AddEntry(gr[1],"Crystalball", "lep");
    l2->AddEntry(gr[5],"gaus", "lp");
    l2->Draw();

	compare[2]->cd();
	compare[2]->SetGrid();
	gr[2]->Draw("AP");
	gr[6]->Draw("P");

	TLegend *l3 = new TLegend(0.6, 0.6, 0.75, 0.75);
    l3->SetHeader("compare3 ","C"); 
    l3->AddEntry(gr[2],"Crystalball", "lep");
    l3->AddEntry(gr[6],"gaus", "lp");
    l3->Draw();


	compare[3]->cd();
	compare[3]->SetGrid();
	gr[3]->Draw("AP");
	gr[7]->Draw("P");

	TLegend *l4 = new TLegend(0.6, 0.6, 0.75, 0.75);
    l4->SetHeader("compare4 ","C"); 
    l4->AddEntry(gr[3],"Crystalball", "lep");
    l4->AddEntry(gr[7],"gaus", "lp");
    l4->Draw();

    // // //compare methods of count peaks//

    // // int num_of_methods = 4;

    // // TH1D *CompareHist [num_of_methods];
    // // TCanvas *CompareCanv[num_of_methods - 1];

    // // for (int k = 0; k < num_of_methods; ++k)
    // // {	
    // // 	CompareHist[k] = new TH1D(Form("CompareHist*%d", k+1), Form("CompareHist*%d", k+1), number_of_proj, 0, number_of_proj);

    // // 	if (k != 4)
    // // 	{
    // // 		CompareCanv[k] = new TCanvas(Form("CompareCanv%d", k+1), Form("CompareCanv%d", k+1));
    // // 	}
    	

    // // 	if (k == 0)
    // // 	{
    // // 		for (int i = 0; i < number_of_proj; ++i)
    // // 		{
    // // 			CompareHist[k] -> SetBinContent(i+1, number_of_pairs[i]);
    // // 		}
    // // 	}
    // // 	if (k == 1)
    // // 	{
    // // 		for (int i = 0; i < number_of_proj; ++i)
    // // 		{
    // // 			CompareHist[k] -> SetBinContent(i+1, number_of_pairs2[i]);
    // // 		}
    // // 	}
    // // 	if (k == 2)
    // // 	{
    // // 		for (int i = 0; i < number_of_proj; ++i)
    // // 		{
    // // 			CompareHist[k] -> SetBinContent(i+1, number_of_pairsGaus[i]);
    // // 		}
    // // 	}
    // // 	if (k == 3)
    // // 	{
    // // 		for (int i = 0; i < number_of_proj; ++i)
    // // 		{
    // // 			CompareHist[k] -> SetBinContent(i+1, number_of_pairs2Gaus[i]);
    // // 		}
    // // 	}
    // // }

    // // for (int i = 1; i < num_of_methods; ++i)
    // // {
    // // 	CompareHist[i]->Divide(CompareHist[0]);

    // // 	CompareCanv[i-1]->cd();
    // // 	CompareCanv[i-1]->SetGrid();

    // // 	CompareHist[i]->Draw();
    // // }

    //double Bins[] = {1.0, 1.5, 2.0, 3.0, 4.0, 4.5, 10.0};
    
    double Bins[] = {1.0, 1.4, 2.0, 2.5, 3.2, 5.0, 10.0};

	TH1D *PeakDataNoField  = new TH1D("PeakDataNoField (No Field)", " ", number_of_proj, Bins);
	TH1D *PeakErDataNoField  = new TH1D("PeakErDataNoField (No Field)", " ", number_of_proj, Bins);

	TH1D *SigmaDataNoField  = new TH1D("SigmaDataNoField (No Field)", " ", number_of_proj, Bins);
	TH1D *SigmaErDataNoField  = new TH1D("SigmaErDataNoField (No Field)", " ", number_of_proj, Bins);

	for (int i = 0; i < number_of_proj; ++i)
	{
		PeakDataNoField ->SetBinContent(i+1, peaks[i]);
		PeakErDataNoField ->SetBinContent(i+1, peaks_er[i]);

		SigmaDataNoField ->SetBinContent(i+1, sigma[i]);
		SigmaErDataNoField ->SetBinContent(i+1, sigma_er[i]);
	}

	// TFile *Compare = new TFile("compare_for_NIRS.root", "update");
	// PeakDataNoField ->Write();
	// PeakErDataNoField ->Write();
	// SigmaDataNoField ->Write();
	// SigmaErDataNoField ->Write();
	// Compare->Close();

    TH1D *dHnormed_number_of_pairsDataNoField  = new TH1D("dnormed_number_of_pairsDataNoField", " ", number_of_proj, Bins);
    TH1D *dHnormed_number_of_pairs2 = new TH1D("CBDataFuncIntegralNoField", " ", number_of_proj,Bins);
    TH1D *dHnormed_number_of_pairsGaus = new TH1D("GausDataHistIntegralNoField", " ", number_of_proj, Bins);
    TH1D *dHnormed_number_of_pairs2Gaus = new TH1D("GausDataFuncIntegralNoField", " ", number_of_proj, Bins);

    TH1D *dHnormed_number_of_pairs_ErDataNoField   = new TH1D("dnormed_number_of_pairs_erDataNoField", " ", number_of_proj, Bins);
    TH1D *dHnormed_number_of_pairs2_er = new TH1D("ErCBDataFuncIntegralNoField", " ", number_of_proj, Bins);
    TH1D *dHnormed_number_of_pairs_erGaus = new TH1D("ErGausDataHistIntegralNoField", " ", number_of_proj, Bins);
    TH1D *dHnormed_number_of_pairs2_erGaus = new TH1D("ErGausDataFuncIntegralNoField", " ", number_of_proj,Bins);

    for (int i = 0; i < number_of_proj; ++i)
    {
    	dHnormed_number_of_pairsDataNoField  ->SetBinContent(i+1, normed_number_of_pairs[i]);
    	dHnormed_number_of_pairs2 ->SetBinContent(i+1, normed_number_of_pairs2[i]);
    	dHnormed_number_of_pairsGaus ->SetBinContent(i+1, normed_number_of_pairsGaus[i]);
    	dHnormed_number_of_pairs2Gaus ->SetBinContent(i+1, normed_number_of_pairs2Gaus[i]);

    	dHnormed_number_of_pairs_ErDataNoField   ->SetBinContent(i+1, sqrt(number_of_pairs_er[i])/normed_koeff);
    	//dHnormed_number_of_pairs_ErDataNoField   ->SetBinContent(i+1, normed_number_of_pairs_er[i]);
    	dHnormed_number_of_pairs2_er ->SetBinContent(i+1, normed_number_of_pairs_er[i]);
    	dHnormed_number_of_pairs_erGaus ->SetBinContent(i+1, normed_number_of_pairs_erGaus[i]);
    	dHnormed_number_of_pairs2_erGaus ->SetBinContent(i+1, normed_number_of_pairs_er2Gaus[i]);

    }

    for (int i = 0; i < number_of_proj; ++i)
    {
    	cout << "\t" <<dHnormed_number_of_pairsDataNoField->GetBinContent(i+1) << "\t" << dHnormed_number_of_pairs_ErDataNoField->GetBinContent(i+1) << "\t" << sqrt(number_of_pairs_er[i])/normed_koeff << endl;
    }

    // TFile *CompareDataMC = new TFile("CompareDataMCNoField.root", "recreate");
    // PeakDataNoField->Write();
    // PeakErDataNoField->Write();
    // SigmaDataNoField->Write();
    // SigmaErDataNoField->Write();
    // dHnormed_number_of_pairsDataNoField->Write();
    // dHnormed_number_of_pairs_ErDataNoField ->Write();
    // CompareDataMC->Close();


    // TFile *ThreeMetods = new TFile("ThreeMetods.root", "update");
    // ThreeMetods->cd("GausHistIntegral");
    // dHnormed_number_of_pairsGaus->Write();
    // dHnormed_number_of_pairs_erGaus->Write();
    // ThreeMetods->Close();

    //TFile *ThreeMetods = new TFile("ThreeMetods.root", "update");
    // ThreeMetods->cd("GausFuncIntegral");
    // dHnormed_number_of_pairs2Gaus->Write();
    // dHnormed_number_of_pairs2_erGaus->Write();
    // ThreeMetods->Close();

    // ThreeMetods->cd("CBFuncIntegral");
    // dHnormed_number_of_pairs2->Write();
    // dHnormed_number_of_pairs2_er->Write();
    // ThreeMetods->Close();

}


		// if (i == 0)
		// {
		// 	func[i]->SetParameters(max_bin[i]*2,
		// 	0.135,
		// 	0.006,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }
		// else if (i == 1)
		// {
		// 	func[i]->SetParameters(max_bin[i]/2,
		// 	0.135,
		// 	0.006,
		// 	90,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);	
		// }

		// else if (i > 1)
		// {
		// 	func[i]->SetParameters(max_bin[i],
		// 	0.135,
		// 	0.006,
		// 	0.05,
		// 	1.7,
		// 	max_bin[i] / 40.,
		// 	-max_bin[i]/ 10.,
		// 	max_bin[i]/ 10.);
		// }

		// cout << "**********************" << endl;
		// cout << i+1  << " -я проекция " << low_pt[i] << " Pt < " << up_pt[i] << endl;
		// cout << "**********************" << endl;

		// real_mixed[i]->Fit("func", "N", "", 0.06, 0.40);

	    // TFile *ResultFile = new TFile("Result_LHC17c.root", "update");
	    // dHnormed_number_of_pairsDataNoField ->Write();
	    // dHnormed_number_of_pairs_ErDataNoField  ->Write();
	    // ResultFile->Close();




	    // TFile *Compare = new TFile("compare_for_NIRS_remake.root", "update");
		// PeakDataNoField ->Write();
		// PeakErDataNoField ->Write();
		// SigmaDataNoField ->Write();
		// SigmaErDataNoField ->Write();
		// Compare->Close();

	    // TFile *ResultFile = new TFile("Result_LHC17c_remake.root", "update");
	    // dHnormed_number_of_pairsDataNoField ->Write();
	    // dHnormed_number_of_pairs_ErDataNoField  ->Write();
	    // ResultFile->Close();
