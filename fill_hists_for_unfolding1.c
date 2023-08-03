#define fill_hists_for_unfolding1_cxx
#include "fill_hists_for_unfolding1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "histio.c"

#include "draw_response.c"

void fill_hists_for_unfolding1::Loop()
{

   bool verbose = false ;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   printf("\n\n Ntuple has %lld entries\n\n", nentries ) ;

   gDirectory -> Delete( "h*" ) ;



 //-------------
// float simple_min = 0.10 ;
// float simple_max = 0.80 ;

// int   simple_nbins_obs = 15 ;
// int   simple_nbins_gen = 10 ;
 //-------------

 //-------------
   float simple_min = 0.10 ;
   float simple_max = 0.80 ;

   int   simple_nbins_obs = 20 ;
   int   simple_nbins_gen = 15 ;
 //-------------

   float adv_min = simple_min ;
   float adv_max = simple_max ;
   int   adv_nbins_obs = simple_nbins_obs ;
   int   adv_nbins_gen = simple_nbins_gen ;

   TH2F* h_simple_gen_vs_obs = new TH2F( "h_simple_gen_vs_obs", "h_simple_gen_vs_obs", simple_nbins_obs, simple_min, simple_max, simple_nbins_gen, simple_min, simple_max ) ;
   TH2F* h_adv_gen_vs_obs = new TH2F( "h_adv_gen_vs_obs", "h_adv_gen_vs_obs", adv_nbins_obs, adv_min, adv_max, adv_nbins_gen, adv_min, adv_max ) ;
   TH2F* h_adv0_gen_vs_obs = new TH2F( "h_adv0_gen_vs_obs", "h_adv0_gen_vs_obs", adv_nbins_obs, adv_min, adv_max, adv_nbins_gen, adv_min, adv_max ) ;



   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      int ei = jentry ;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if ( !verbose && ei%100 == 0 ) {
         printf(" --- Event: %7d / %lld    %6.3f\r", ei, nentries, (1.*ei)/(1.*nentries) ) ;
         fflush(stdout) ;
      }

      h_simple_gen_vs_obs -> Fill( simple_classifier_output_obs, simple_classifier_output_gen ) ;
      h_adv_gen_vs_obs -> Fill( classifier_output_obs, classifier_output_gen ) ;
      h_adv0_gen_vs_obs -> Fill( classifier_output_obs_lam0, classifier_output_gen_lam0 ) ;


   } // jentry


   char outfile[1000] ;
   sprintf( outfile, "unfold-hists-%s", filename ) ;
   saveHist( outfile, "h*" ) ;

   draw_response() ;

}
