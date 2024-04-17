

//


#include "histio.c"
#include "utils.c"

#include "tunfold-17.9/TUnfoldDensity.h"

//----------

   void chi2_study1a(
                 const char* hist_name = "h_sf3_gen_vs_obs",
                 int rn_seed = 1234,
                 int ntoy = 100,
                 int ngen = 2000,
                 const char* input_fileA = "unfold-hists2-binning-24-12-output_bce_mse_rapgap.root",
                 const char* input_fileB = "unfold-hists2-binning-24-12-output_bce_mse_django.root"
                 ) {

      bool verbose(false) ;

      gRandom -> SetSeed( rn_seed ) ;

      printf("\n\n Loading TUnfold shared library.\n\n") ;
      gSystem -> Load( "tunfold-17.9/libunfold.a" ) ;
      printf("\n\n Done.\n\n") ;

      gDirectory -> Delete( "*" ) ;


      loadHist( input_fileA, "A" ) ;
      loadHist( input_fileB, "B" ) ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      TH2F* h_in_gen_vs_obsA = get_hist2d( Form("%s_A", hist_name) ) ;
      TH2F* h_in_gen_vs_obsB = get_hist2d( Form("%s_B", hist_name) ) ;


      TH1* h_obs_sourceA = h_in_gen_vs_obsA -> ProjectionX( "h_obs_sourceA" ) ;
      TH1* h_gen_sourceA = h_in_gen_vs_obsA -> ProjectionY( "h_gen_sourceA" ) ;

      TH1* h_obs_sourceB = h_in_gen_vs_obsB -> ProjectionX( "h_obs_sourceB" ) ;
      TH1* h_gen_sourceB = h_in_gen_vs_obsB -> ProjectionY( "h_gen_sourceB" ) ;


      int nbins = 40 ;
      double hchi2max = 250. ;
      TH1F* h_chi2_ufA = new TH1F( "h_chi2_ufA", "chi2 unfolding only A", nbins, 0., hchi2max ) ;
      TH1F* h_chi2_ufB = new TH1F( "h_chi2_ufB", "chi2 unfolding only B", nbins, 0., hchi2max ) ;
      TH1F* h_chi2_totalA = new TH1F( "h_chi2_totalA", "chi2 total A", nbins, 0., hchi2max ) ;
      TH1F* h_chi2_totalB = new TH1F( "h_chi2_totalB", "chi2 total B", nbins, 0., hchi2max ) ;


      for ( int ti=0; ti<ntoy; ti++ ) {

         gDirectory -> Delete( "ht_*" ) ;

         TH1* h_obs_randomA = (TH1*) h_obs_sourceA -> Clone( "ht_obs_randomA" ) ;
         h_obs_randomA->Reset() ;
         h_obs_randomA->FillRandom( h_obs_sourceA, ngen ) ;

         TH1* h_obs_randomB = (TH1*) h_obs_randomA -> Clone( "ht_obs_randomB" ) ;
         /////////h_obs_randomB->Reset() ;
         /////////h_obs_randomB->FillRandom( h_obs_sourceB, ngen ) ;


         if (verbose) printf("\n\n ----- Creating a TUnfoldDensity instance now for A.\n") ;
         TUnfoldDensity unfoldA( h_in_gen_vs_obsA, TUnfold::kHistMapOutputVert ) ;
         if (verbose) printf("\n\n ----- Done.\n") ;

         if (verbose) printf("\n\n ----- Creating a TUnfoldDensity instance now for B.\n") ;
         TUnfoldDensity unfoldB( h_in_gen_vs_obsB, TUnfold::kHistMapOutputVert ) ;
         if (verbose) printf("\n\n ----- Done.\n") ;


         int return_statusA = unfoldA.SetInput( h_obs_randomA ) ;
         if (verbose) printf("  Return status for SetInput A: %d\n", return_statusA ) ;

         int return_statusB = unfoldB.SetInput( h_obs_randomB ) ;
         if (verbose) printf("  Return status for SetInput B: %d\n", return_statusB ) ;



         //========================================================================
         // the unfolding is done here
         //
         // scan L curve and find best point
         Int_t nScan=30;
         // use automatic L-curve scan: start with taumin=taumax=0.0
         Double_t tauMin=0.0;
         Double_t tauMax=0.0;

         Int_t iBestA, iBestB;
         TSpline *logTauXA, *logTauXB, *logTauYA, *logTauYB;
         TGraph *lCurveA, *lCurveB;

         // if required, report Info messages (for debugging the L-curve scan)
         Int_t oldinfo=gErrorIgnoreLevel;
         gErrorIgnoreLevel=kInfo;

         // this method scans the parameter tau and finds the kink in the L curve
         // finally, the unfolding is done for the best choice of tau
         if (verbose) printf("\n Begin unfolding for A.\n") ;
         iBestA=unfoldA.ScanLcurve(nScan,tauMin,tauMax,&lCurveA,&logTauXA,&logTauYA);
         if (verbose) printf("\n Begin unfolding for B.\n") ;
         iBestB=unfoldB.ScanLcurve(nScan,tauMin,tauMax,&lCurveB,&logTauXB,&logTauYB);

         // if required, switch to previous log-level
         gErrorIgnoreLevel=oldinfo;


         //==========================================================================
         // print some results
         //
         if (verbose) {
            printf("\n Results for A.\n") ;
            std::cout<<"tau="<<unfoldA.GetTau()<<"\n";
            std::cout<<"chi**2="<<unfoldA.GetChi2A()<<"+"<<unfoldA.GetChi2L()
                  <<" / "<<unfoldA.GetNdf()<<"\n";
         }




         if (verbose) {
            printf("\n Results for B.\n") ;
            std::cout<<"tau="<<unfoldB.GetTau()<<"\n";
            std::cout<<"chi**2="<<unfoldB.GetChi2A()<<"+"<<unfoldB.GetChi2L()
                  <<" / "<<unfoldB.GetNdf()<<"\n";
         }


         //==========================================================================
         // retrieve results into histograms

         // get unfolded distribution
         TH1 *histMunfoldA=unfoldA.GetOutput("ht_UnfoldedA");
         TH1 *histMunfoldB=unfoldB.GetOutput("ht_UnfoldedB");

         // get unfolding result, folded back
         TH1 *histMdetFoldA=unfoldA.GetFoldedOutput("ht_FoldedBackA");
         TH1 *histMdetFoldB=unfoldB.GetFoldedOutput("ht_FoldedBackB");

         // get error matrix (input distribution [stat] errors only)
         // TH2D *histEmatData=unfold.GetEmatrix("EmatData");

         // get total error matrix:
         //   migration matrix uncorrelated and correlated systematic errors
         //   added in quadrature to the data statistical errors
         TH2 *histEmatTotalA=unfoldA.GetEmatrixTotal("ht_EmatTotalA");
         TH2 *histEmatTotalB=unfoldB.GetEmatrixTotal("ht_EmatTotalB");

         // create data histogram with the total errors
         int nGen = histMunfoldA -> GetNbinsX() ;
         float xminGen = histMunfoldA -> GetXaxis() -> GetXmin() ;
         float xmaxGen = histMunfoldA -> GetXaxis() -> GetXmax() ;

         TH1D *histTotalErrorA= new TH1D("ht_TotalErrorA","TotalErrorA",nGen,xminGen,xmaxGen);
         TH1D *histTotalErrorB= new TH1D("ht_TotalErrorB","TotalErrorB",nGen,xminGen,xmaxGen);

         for(Int_t bin=1;bin<=nGen;bin++) {
           histTotalErrorA->SetBinContent(bin,histMunfoldA->GetBinContent(bin));
           histTotalErrorB->SetBinContent(bin,histMunfoldB->GetBinContent(bin));
           histTotalErrorA->SetBinError (bin,TMath::Sqrt(histEmatTotalA->GetBinContent(bin,bin)));
           histTotalErrorB->SetBinError (bin,TMath::Sqrt(histEmatTotalB->GetBinContent(bin,bin)));
         }



         //-- calculate the correlation coefficients matrix
         TH2* correlation_matrixA = (TH2*) histEmatTotalA -> Clone( "ht_correlation_matrixA" ) ;
         TH2* correlation_matrixB = (TH2*) histEmatTotalB -> Clone( "ht_correlation_matrixB" ) ;

         correlation_matrixA -> SetTitle( "Correlation coefficients A" ) ;
         correlation_matrixB -> SetTitle( "Correlation coefficients B" ) ;

         for ( int xbi=1; xbi<=histEmatTotalA->GetNbinsX(); xbi++ ) {
            for ( int ybi=1; ybi<=histEmatTotalA->GetNbinsY(); ybi++ ) {

               float rho = 1. ;
               float sigma_x2, sigma_y2 ;

               sigma_x2 = histEmatTotalA->GetBinContent( xbi, xbi ) ;
               sigma_y2 = histEmatTotalA->GetBinContent( ybi, ybi ) ;
               if ( sigma_x2 > 0 && sigma_y2 > 0 ) {
                  rho = histEmatTotalA -> GetBinContent( xbi, ybi )  / sqrt( sigma_x2 * sigma_y2 ) ;
               }
               correlation_matrixA -> SetBinContent( xbi, ybi, rho ) ;

               sigma_x2 = histEmatTotalB->GetBinContent( xbi, xbi ) ;
               sigma_y2 = histEmatTotalB->GetBinContent( ybi, ybi ) ;
               if ( sigma_x2 > 0 && sigma_y2 > 0 ) {
                  rho = histEmatTotalB -> GetBinContent( xbi, ybi )  / sqrt( sigma_x2 * sigma_y2 ) ;
               }
               correlation_matrixB -> SetBinContent( xbi, ybi, rho ) ;

            } // ybi
         } // xbi

         correlation_matrixA -> SetMinimum( -1. ) ;
         correlation_matrixB -> SetMinimum( -1. ) ;


         // get global correlation coefficients
         // for this calculation one has to specify whether the
         // underflow/overflow bins are included or not
         // default: include all bins
         // here: exclude underflow and overflow bins
         TH2 *gHistInvEMatrixA, *gHistInvEMatrixB ;
         TH1 *histRhoiA=unfoldA.GetRhoItotal("ht_rho_IA",
                                           0, // use default title
                                           0, // all distributions
                                           "*[UO]", // discard underflow and overflow bins on all axes
                                           kTRUE, // use original binning
                                           &gHistInvEMatrixA // store inverse of error matrix
                                           );
         TH1 *histRhoiB=unfoldB.GetRhoItotal("ht_rho_IB",
                                           0, // use default title
                                           0, // all distributions
                                           "*[UO]", // discard underflow and overflow bins on all axes
                                           kTRUE, // use original binning
                                           &gHistInvEMatrixB // store inverse of error matrix
                                           );


         TH2F* h_normalized_responseA = (TH2F*) h_in_gen_vs_obsA -> Clone( "ht_normalized_responseA" ) ;
         TH2F* h_normalized_responseB = (TH2F*) h_in_gen_vs_obsB -> Clone( "ht_normalized_responseB" ) ;

         for ( int ybi=1; ybi<=h_in_gen_vs_obsA->GetNbinsY(); ybi++ ) {
            float row_sumA = 0. ;
            float row_sumB = 0. ;
            for ( int xbi=1; xbi<=h_in_gen_vs_obsA->GetNbinsX(); xbi++ ) {
               row_sumA += h_in_gen_vs_obsA -> GetBinContent( xbi, ybi ) ;
               row_sumB += h_in_gen_vs_obsB -> GetBinContent( xbi, ybi ) ;
            } // xbi
            for ( int xbi=1; xbi<=h_in_gen_vs_obsA->GetNbinsX(); xbi++ ) {
               if ( row_sumA > 0 ) {
                  h_normalized_responseA -> SetBinContent( xbi, ybi, ( h_in_gen_vs_obsA -> GetBinContent( xbi, ybi ) ) / row_sumA ) ;
                  h_normalized_responseA -> SetBinError( xbi, ybi, 0. ) ;
               } else {
                  h_normalized_responseA -> SetBinContent( xbi, ybi, 0. ) ;
                  h_normalized_responseA -> SetBinError( xbi, ybi, 0. ) ;
               }
               if ( row_sumB > 0 ) {
                  h_normalized_responseB -> SetBinContent( xbi, ybi, ( h_in_gen_vs_obsB -> GetBinContent( xbi, ybi ) ) / row_sumB ) ;
                  h_normalized_responseB -> SetBinError( xbi, ybi, 0. ) ;
               } else {
                  h_normalized_responseB -> SetBinContent( xbi, ybi, 0. ) ;
                  h_normalized_responseB -> SetBinError( xbi, ybi, 0. ) ;
               }
            } // xbi
         } // ybi


         TH1* h_gen_compareA = (TH1*) h_gen_sourceA -> Clone( "ht_gen_compareA" ) ;
         TH1* h_gen_compareB = (TH1*) h_gen_sourceB -> Clone( "ht_gen_compareB" ) ;
         h_gen_compareA -> Scale( ( histMunfoldA -> Integral() )/( h_gen_sourceA -> Integral() ) ) ;
         h_gen_compareB -> Scale( ( histMunfoldB -> Integral() )/( h_gen_sourceB -> Integral() ) ) ;



       //---- Compute vector of differences
         TH1* h_diff_uf_plA = (TH1*) histMunfoldA -> Clone( "ht_diff_uf_plA" ) ;
         ////////printf("just before add a\n") ;
         ////////h_gen_compareA->Print("all") ;
         ////////h_diff_uf_plA->Print("all") ;
         h_diff_uf_plA -> Add( h_gen_compareA, -1. ) ;
         ////////printf("just after add a\n") ;
         TH1* h_diff_uf_plB = (TH1*) histMunfoldB -> Clone( "ht_diff_uf_plB" ) ;
         h_diff_uf_plB -> Add( h_gen_compareB, -1. ) ;

         TVectorD vec_diff_uf_plA( h_diff_uf_plA -> GetNbinsX() ) ;
         TVectorD vec_diff_uf_plB( h_diff_uf_plB -> GetNbinsX() ) ;
         for ( int bi=1; bi<=h_diff_uf_plA -> GetNbinsX(); bi++ ) {
            vec_diff_uf_plA(bi-1) = h_diff_uf_plA -> GetBinContent(bi) ;
            vec_diff_uf_plB(bi-1) = h_diff_uf_plB -> GetBinContent(bi) ;
         }


       //---- Fill covariance matrix from output histogram
         TMatrixD mat_cov_ufA( h_diff_uf_plA -> GetNbinsX(), h_diff_uf_plA -> GetNbinsX() ) ;
         TMatrixD mat_cov_ufB( h_diff_uf_plB -> GetNbinsX(), h_diff_uf_plB -> GetNbinsX() ) ;

         for ( int xbi=1; xbi<=h_diff_uf_plA -> GetNbinsX(); xbi++ ) {
            for ( int ybi=1; ybi<=h_diff_uf_plA -> GetNbinsX(); ybi++ ) {
               mat_cov_ufA( xbi-1, ybi-1 ) = histEmatTotalA -> GetBinContent( xbi, ybi ) ;
               mat_cov_ufB( xbi-1, ybi-1 ) = histEmatTotalB -> GetBinContent( xbi, ybi ) ;
            } // ybi
         } // xbi



       //---- Invert covariance matrix
         TMatrixD mat_inv_cov_ufA( h_diff_uf_plA -> GetNbinsX(), h_diff_uf_plA -> GetNbinsX() ) ;
         mat_inv_cov_ufA = mat_cov_ufA ;
         mat_inv_cov_ufA.Invert() ;

         TMatrixD mat_inv_cov_ufB( h_diff_uf_plB -> GetNbinsX(), h_diff_uf_plB -> GetNbinsX() ) ;
         mat_inv_cov_ufB = mat_cov_ufB ;
         mat_inv_cov_ufB.Invert() ;

         if (verbose) printf("\n\n") ;
         if (verbose) printf(" Output histogram for covariance matrix.\n") ;
         if (verbose) histEmatTotalA->Print("all") ;
         if (verbose) printf("\n\n") ;
         if (verbose) printf(" Covariance matrix A\n") ;
         if (verbose) mat_cov_ufA.Print() ;
         if (verbose) printf("\n\n") ;
         if (verbose) printf(" Inverted Covariance matrix A\n") ;
         if (verbose) mat_inv_cov_ufA.Print() ;




        //---- Compute chi2

         TVectorT<double> prod1 = mat_inv_cov_ufA * vec_diff_uf_plA ;
         if (verbose) printf(" Inverted cov mat times diff vec:\n") ;
         if (verbose) prod1.Print() ;

         double prod2 = ( vec_diff_uf_plA ) * prod1 ;
         if (verbose) printf("  diff vec times prod1 :  %.9f\n", prod2 ) ;

         double chi2A = ( vec_diff_uf_plA ) * ( mat_inv_cov_ufA * vec_diff_uf_plA ) ;
         if (verbose) printf("  chi2 A : %.9f\n", chi2A ) ;

         double chi2B = ( vec_diff_uf_plB ) * ( mat_inv_cov_ufB * vec_diff_uf_plB ) ;
         if (verbose) printf("  chi2 B : %.9f\n", chi2B ) ;






        //---- Calculate the covariance matrix from the unfolding model systematic.

         TMatrixD mat_cov_model_syst( h_diff_uf_plA -> GetNbinsX(), h_diff_uf_plA -> GetNbinsX() ) ;

         for ( int xbi=1; xbi<=h_diff_uf_plA -> GetNbinsX(); xbi++ ) {
            for ( int ybi=1; ybi<=h_diff_uf_plA -> GetNbinsX(); ybi++ ) {
               mat_cov_model_syst( xbi-1, ybi-1 ) = ( histMunfoldB->GetBinContent(xbi) - histMunfoldA->GetBinContent(xbi) ) *
                                                    ( histMunfoldB->GetBinContent(ybi) - histMunfoldA->GetBinContent(ybi) ) ;
            } // ybi
         } // xbi

         if (verbose) printf("\n\n") ;
         if (verbose) printf("  Covariance matrix for unfolding model systematic.\n") ;
         if (verbose) mat_cov_model_syst.Print() ;




        //--- combine unfolding and model syst covariance matrices
         TMatrixD mat_cov_totalA = mat_cov_ufA + mat_cov_model_syst ;
         TMatrixD mat_cov_totalB = mat_cov_ufB + mat_cov_model_syst ;

         if (verbose) printf("  Total covariance matrix A: \n" ) ;
         if (verbose) mat_cov_totalA.Print() ;

         TMatrixD mat_inv_cov_totalA( h_diff_uf_plA -> GetNbinsX(), h_diff_uf_plA -> GetNbinsX() ) ;
         mat_inv_cov_totalA = mat_cov_totalA ;
         mat_inv_cov_totalA.Invert() ;

         TMatrixD mat_inv_cov_totalB( h_diff_uf_plA -> GetNbinsX(), h_diff_uf_plA -> GetNbinsX() ) ;
         mat_inv_cov_totalB = mat_cov_totalB ;
         mat_inv_cov_totalB.Invert() ;

         if (verbose) printf("  Inverted Total covariance matrix A: \n" ) ;
         if (verbose) mat_inv_cov_totalA.Print() ;


         double chi2_totalA = ( vec_diff_uf_plA ) * ( mat_inv_cov_totalA * vec_diff_uf_plA ) ;
         double chi2_totalB = ( vec_diff_uf_plB ) * ( mat_inv_cov_totalB * vec_diff_uf_plB ) ;

         printf("  chi2 A : %.9f,  total chi2 A :  %.9f\n", chi2A, chi2_totalA ) ;
         printf("  chi2 B : %.9f,  total chi2 B :  %.9f\n", chi2B, chi2_totalB ) ;


         h_chi2_ufA -> Fill( chi2A ) ;
         h_chi2_ufB -> Fill( chi2B ) ;

         h_chi2_totalA -> Fill( chi2_totalA ) ;
         h_chi2_totalB -> Fill( chi2_totalB ) ;


      } // ti





      h_chi2_ufA -> SetLineColor(4) ;
      h_chi2_totalA -> SetLineColor(4) ;

      h_chi2_ufB -> SetLineColor(2) ;
      h_chi2_totalB -> SetLineColor(2) ;

      h_chi2_ufA -> SetLineStyle(2) ;
      h_chi2_ufB -> SetLineStyle(2) ;

      h_chi2_ufA -> SetLineWidth(3) ;
      h_chi2_ufB -> SetLineWidth(3) ;
      h_chi2_totalA -> SetLineWidth(3) ;
      h_chi2_totalB -> SetLineWidth(3) ;

      h_chi2_totalA -> SetXTitle( "#chi^{2}" ) ;
      h_chi2_totalA -> SetYTitle( "Toy experiments" ) ;

      h_chi2_totalA -> SetTitleSize( 0.050, "x" ) ;
      h_chi2_totalA -> SetTitleSize( 0.050, "y" ) ;

      h_chi2_totalB -> SetXTitle( "#chi^{2}" ) ;
      h_chi2_totalB -> SetYTitle( "Toy experiments" ) ;

      h_chi2_totalB -> SetTitleSize( 0.050, "x" ) ;
      h_chi2_totalB -> SetTitleSize( 0.050, "y" ) ;

      h_chi2_totalA -> SetTitle( "A compared to A" ) ;
      h_chi2_totalB -> SetTitle( "A compared to B" ) ;


      if ( h_chi2_totalB -> GetMaximumBin() > 5 ) {
         printf("\n Resetting max for h_chi2_totalB\n") ;
         h_chi2_totalB -> SetMaximum( 1.4 * (h_chi2_totalB->GetMaximum()) ) ;
      }


      float lx = 0.38 ;
      float ly = 0.73 ;
      float lw = 0.35 ;
      float lh = 0.22 ;
      TLegend* la = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      TLegend* lb = new TLegend( lx, ly, lx+lw, ly+lh ) ;

      la -> AddEntry( h_chi2_ufA, "Stat. only" ) ;
      la -> AddEntry( h_chi2_totalA, "Stat. + model syst." ) ;

      lb -> AddEntry( h_chi2_ufB, "Stat. only" ) ;
      lb -> AddEntry( h_chi2_totalB, "Stat. + model syst." ) ;

      gStyle -> SetOptStat( "mr" ) ;


      float sx = 0.75 ;
      float sy = 0.85 ;
      float sw = 0.20 ;
      float sh = 0.10 ;
      float sdy = 0.12 ;

      TCanvas* can1 = get_canvas( "can1", "", 50, 50, 800, 1200 ) ;

      can1 -> Clear() ;
      can1 -> Divide(1,2) ;

      can1 -> cd(1) ;
      h_chi2_totalA -> Draw() ;
      can1 -> Update() ;
      TPaveStats* sta1 = (TPaveStats*) h_chi2_totalA -> GetListOfFunctions() -> FindObject("stats") ;
      sta1->SetX1NDC(sx) ;
      sta1->SetX2NDC(sx+sw) ;
      sta1->SetY1NDC(sy-sdy) ;
      sta1->SetY2NDC(sy+sh-sdy) ;
      h_chi2_ufA -> Draw() ;
      can1 -> Update() ;
      TPaveStats* sta2 = (TPaveStats*) h_chi2_ufA -> GetListOfFunctions() -> FindObject("stats") ;
      sta2->SetX1NDC(sx) ;
      sta2->SetX2NDC(sx+sw) ;
      sta2->SetY1NDC(sy) ;
      sta2->SetY2NDC(sy+sh) ;
      h_chi2_totalA->Draw() ;
      h_chi2_ufA->Draw("same") ;
      sta1->Draw() ;
      sta2->Draw() ;
      la -> Draw() ;

      can1 -> cd(2) ;
      h_chi2_totalB -> Draw("same") ;
      can1 -> Update() ;
      TPaveStats* stb1 = (TPaveStats*) h_chi2_totalB -> GetListOfFunctions() -> FindObject("stats") ;
      stb1->SetX1NDC(sx) ;
      stb1->SetX2NDC(sx+sw) ;
      stb1->SetY1NDC(sy-sdy) ;
      stb1->SetY2NDC(sy+sh-sdy) ;
      h_chi2_ufB -> Draw() ;
      can1 -> Update() ;
      TPaveStats* stb2 = (TPaveStats*) h_chi2_ufB -> GetListOfFunctions() -> FindObject("stats") ;
      stb2->SetX1NDC(sx) ;
      stb2->SetX2NDC(sx+sw) ;
      stb2->SetY1NDC(sy) ;
      stb2->SetY2NDC(sy+sh) ;
      h_chi2_totalB->Draw() ;
      h_chi2_ufB->Draw("same") ;
      stb1->Draw() ;
      stb2->Draw() ;
      lb -> Draw() ;






     //===============================================================

///   gStyle -> SetOptStat(0) ;

///   h_obs_sourceA -> SetLineColor(4) ;
///   h_gen_sourceA -> SetLineColor(2) ;

///   h_obs_sourceB -> SetLineColor(4) ;
///   h_gen_sourceB -> SetLineColor(2) ;


///   TExec* change_hist_palette = new TExec( "change_hist_palette", "Setup2DhistPalette();" );
///   TExec* change_cor_palette = new TExec( "change_cor_palette", "SetupCorrelationPalette();" );



///   TCanvas* can1A = (TCanvas*) gDirectory -> FindObject( "can1A" ) ;
///   if ( can1A == 0x0 ) can1A = new TCanvas( "can1A", "A", 50, 50, 1000, 1000 ) ;
///   can1A -> Clear() ;
///   can1A -> cd() ;

///   can1A -> Divide(2,2) ;


///   can1A -> cd(1) ;
///   h_normalized_responseA -> Draw("colz") ;
///   change_hist_palette -> Draw() ;
///   h_normalized_responseA -> Draw("colz same") ;

///   can1A -> cd(2) ;
///   histMunfoldA -> Draw() ;
///   h_gen_compareA -> Draw("same hist") ;

///   can1A -> cd(3) ;

///   correlation_matrixA->Draw("colz") ;
///   change_cor_palette -> Draw() ;
///   correlation_matrixA->Draw("colz same") ;


///   can1A -> cd(4) ;
///   h_gen_sourceA -> Draw("hist") ;








///   TCanvas* can1B = (TCanvas*) gDirectory -> FindObject( "can1B" ) ;
///   if ( can1B == 0x0 ) can1B = new TCanvas( "can1B", "B", 1050, 50, 1000, 1000 ) ;
///   can1B -> Clear() ;
///   can1B -> cd() ;

///   can1B -> Divide(2,2) ;


///   can1B -> cd(1) ;
///   h_normalized_responseB -> Draw("colz") ;
///   change_hist_palette -> Draw() ;
///   h_normalized_responseB -> Draw("colz same") ;

///   can1B -> cd(2) ;
///   histMunfoldB -> Draw() ;
///   h_gen_compareB -> Draw("same hist") ;

///   can1B -> cd(3) ;

///   correlation_matrixB->Draw("colz") ;
///   change_cor_palette -> Draw() ;
///   correlation_matrixB->Draw("colz same") ;


///   can1B -> cd(4) ;
///   h_gen_sourceB -> Draw("hist") ;











   }













