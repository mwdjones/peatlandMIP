CDF      
      lndgrid       gridcell      landunit      column         pft    @   levgrnd       levurb        levlak     
   numrad        levsno        ltype      	   natpft        string_length         levdcmp       levtrc     
   hist_interval         time             title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 11/27/23 16:13:01   source        Community Land Model CLM4.0    hostname      ubuntu     username      mwjones    version       5823c39    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       -TESTShi2015_US-SPR_ICB1850CNRDCTCBC_ad_spinup      Surface_dataset       surfdata.nc    Initial_conditions_dataset        arbitrary initialization   #PFT_physiological_constants_dataset       clm_params.nc      ltype_vegetated_or_bare_soil            
ltype_crop              ltype_landice               (ltype_landice_multiple_elevation_classes            ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   Time_constant_3Dvars_filename         K./TESTShi2015_US-SPR_ICB1850CNRDCTCBC_ad_spinup.clm2.h0.0001-01-01-00000.nc    Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE      !   levgrnd                	long_name         coordinate soil levels     units         m         <     ��   levlak                 	long_name         coordinate lake levels     units         m         (     ��   levdcmp                	long_name         coordinate soil levels     units         m         <     �   time               	long_name         time   units         days since 0001-01-01 00:00:00     calendar      noleap     bounds        time_bounds            �x   mcdate                 	long_name         current date (YYYYMMDD)            �|   mcsec                  	long_name         current seconds of current date    units         s              �   mdcur                  	long_name         current day (from base day)            �   mscur                  	long_name         current seconds of current day             �   nstep                  	long_name         	time step              �   time_bounds                   	long_name         history time interval endpoints            �   date_written                            �   time_written                            �   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��           �@   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��           �H   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��           �P   topo                	long_name         grid cell topography   units         m      
_FillValue        {@��   missing_value         {@��           �X   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��           �`   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����           �h   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����           �p   ACTUAL_IMMOB                   	long_name         actual N immobilization    units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   ACTUAL_IMMOB_P                     	long_name         actual P immobilization    units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   ADSORBTION_P                   	long_name         adsorb P flux      units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   AGNPP                      	long_name         aboveground NPP    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   AGWDNPP                    	long_name         aboveground wood NPP   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   ALT                    	long_name         current active layer thickness     units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   ALTMAX                     	long_name         %maximum annual active layer thickness      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   ALTMAX_LASTYEAR                    	long_name         )maximum prior year active layer thickness      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   AR                     	long_name         !autotrophic respiration (MR + GR)      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    AVAILC                     	long_name         C flux available for allocation    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   AVAIL_RETRANSP                     	long_name         *P flux available from retranslocation pool     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   BAF_CROP                   	long_name         fractional area burned for crop    units         proportion/sec     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	BAF_PEATF                      	long_name         "fractional area burned in peatland     units         proportion/sec     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    BCDEP                      	long_name         -total BC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �(   BGNPP                      	long_name         belowground NPP    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �0   BIOCHEM_PMIN                   	long_name         $biochemical rate of P mineralization   units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �8   BTRAN                      	long_name         transpiration beta factor      units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   	BUILDHEAT                      	long_name         8heat flux from urban building interior to walls and roof   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   CH4PROD                    	long_name          Gridcell total production of CH4   units         gC/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �P   CH4_SURF_AERE_SAT                      	long_name         :aerenchyma surface CH4 flux for inundated area; (+ to atm)     units         mol/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �X   CH4_SURF_AERE_UNSAT                    	long_name         >aerenchyma surface CH4 flux for non-inundated area; (+ to atm)     units         mol/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �`   CH4_SURF_DIFF_SAT                      	long_name         @diffusive surface CH4 flux for inundated / lake area; (+ to atm)   units         mol/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �h   CH4_SURF_DIFF_UNSAT                    	long_name         =diffusive surface CH4 flux for non-inundated area; (+ to atm)      units         mol/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �p   CH4_SURF_EBUL_SAT                      	long_name         Aebullition surface CH4 flux for inundated / lake area; (+ to atm)      units         mol/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �x   CH4_SURF_EBUL_UNSAT                    	long_name         >ebullition surface CH4 flux for non-inundated area; (+ to atm)     units         mol/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   
COL_PTRUNC                     	long_name         "column-level sink for P truncation     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CONC_CH4_SAT                      	long_name         0CH4 soil Concentration for inundated / lake area   units         mol/m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   CONC_CH4_UNSAT                        	long_name         -CH4 soil Concentration for non-inundated area      units         mol/m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   CONC_O2_SAT                       	long_name         /O2 soil Concentration for inundated / lake area    units         mol/m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   CONC_O2_UNSAT                         	long_name         ,O2 soil Concentration for non-inundated area   units         mol/m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   CPOOL                      	long_name         temporary photosynthate C pool     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �p   CWDC                   	long_name         CWD C      units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �x   CWDC_HR                    	long_name         /coarse woody debris C heterotrophic respiration    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	CWDC_LOSS                      	long_name         coarse woody debris C loss     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CWDC_TO_LITR2C                     	long_name         .decomp. of coarse woody debris C to litter 2 C     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CWDC_TO_LITR3C                     	long_name         .decomp. of coarse woody debris C to litter 3 C     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CWDC_vr                       	long_name         CWD C (vertically resolved)    units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   CWDN                   	long_name         CWD N      units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CWDN_TO_LITR2N                     	long_name         .decomp. of coarse woody debris N to litter 2 N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    CWDN_TO_LITR3N                     	long_name         .decomp. of coarse woody debris N to litter 3 N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �(   CWDN_vr                       	long_name         CWD N (vertically resolved)    units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �0   CWDP                   	long_name         CWD P      units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CWDP_TO_LITR2P                     	long_name         .decomp. of coarse woody debris P to litter 2 N     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CWDP_TO_LITR3P                     	long_name         .decomp. of coarse woody debris P to litter 3 N     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   CWDP_vr                       	long_name         CWD P (vertically resolved)    units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   
DEADCROOTC                     	long_name         dead coarse root C     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �8   
DEADCROOTN                     	long_name         dead coarse root N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   
DEADCROOTP                     	long_name         dead coarse root P     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   	DEADSTEMC                      	long_name         dead stem C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �P   	DEADSTEMN                      	long_name         dead stem N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �X   	DEADSTEMP                      	long_name         dead stem P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �`   DENIT                      	long_name         total rate of denitrification      units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �h   DESORPTION_P                   	long_name         desorp P flux      units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �p   DISPVEGC                   	long_name         1displayed veg carbon, excluding storage and cpool      units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �x   DISPVEGN                   	long_name         displayed vegetation nitrogen      units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DISPVEGP                   	long_name         displayed vegetation phosphorus    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DSTDEP                     	long_name         /total dust deposition (dry+wet) from atmosphere    units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DSTFLXT                    	long_name         total surface dust emission    units         kg/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DWB                    	long_name         net change in total water mass     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DWT_CONV_CFLUX_DRIBBLED                    	long_name         Gconversion C flux (immediate loss to atm), dribbled throughout the year    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DWT_CONV_CFLUX_GRC                     	long_name         Xconversion C flux (immediate loss to atm) (0 at all times except first timestep of year)   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DWT_CONV_NFLUX_GRC                     	long_name         Xconversion C flux (immediate loss to atm) (0 at all times except first timestep of year)   units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   DWT_CONV_PFLUX_GRC                     	long_name         Xconversion C flux (immediate loss to atm) (0 at all times except first timestep of year)   units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   DWT_SLASH_CFLUX                    	long_name         .slash C flux to litter and CWD due to land use     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   DWT_SLASH_NFLUX                    	long_name         .slash N flux to litter and CWD due to land use     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   DWT_SLASH_PFLUX                    	long_name         .slash P flux to litter and CWD due to land use     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   EFLX_DYNBAL                    	long_name         0dynamic land cover change conversion energy flux   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   EFLX_GRND_LAKE                     	long_name         Bnet heat flux into lake/snow surface, excluding light transmission     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   EFLX_LH_TOT                    	long_name         !total latent heat flux [+ to atm]      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   EFLX_LH_TOT_R                      	long_name         Rural total evaporation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   EFLX_LH_TOT_U                      	long_name         Urban total evaporation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    ELAI                   	long_name         !exposed one-sided leaf area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   ER                     	long_name         8total ecosystem respiration, autotrophic + heterotrophic   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   ERRH2O                     	long_name         total water conservation error     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	ERRH2OSNO                      	long_name         &imbalance in snow depth (liquid water)     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    ERRSEB                     	long_name         !surface energy conservation error      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �(   ERRSOI                     	long_name         #soil/lake energy conservation error    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �0   ERRSOL                     	long_name         "solar radiation conservation error     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �8   ESAI                   	long_name         !exposed one-sided stem area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   FAREA_BURNED                   	long_name         timestep fractional area burned    units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   FCEV                   	long_name         canopy evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �P   FCH4                   	long_name         2Gridcell surface CH4 flux to atmosphere (+ to atm)     units         kgC/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �X   	FCH4TOCO2                      	long_name          Gridcell oxidation of CH4 to CO2   units         gC/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �`   
FCH4_DFSAT                     	long_name         BCH4 additional flux due to changing fsat, vegetated landunits only     units         kgC/m2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �h   FCOV                   	long_name         fractional impermeable area    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �p   FCTR                   	long_name         canopy transpiration   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �x   FGEV                   	long_name         ground evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FGR                    	long_name         Oheat flux into soil/snow including snow melt and lake / snow light transmission    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FGR12                      	long_name         %heat flux between soil layers 1 and 2      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FGR_R                      	long_name         NRural heat flux into soil/snow including snow melt and snow light transmission     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FGR_U                      	long_name         2Urban heat flux into soil/snow including snow melt     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FH2OSFC                    	long_name         +fraction of ground covered by surface water    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   
FINUNDATED                     	long_name         .fractional inundated area of vegetated columns     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FINUNDATED_LAG                     	long_name         3time-lagged inundated fraction of vegetated columns    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FIRA                   	long_name         !net infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FIRA_R                     	long_name         'Rural net infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FIRA_U                     	long_name         'Urban net infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FIRE                   	long_name         %emitted infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FIRE_R                     	long_name         +Rural emitted infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FIRE_U                     	long_name         +Urban emitted infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FLDS                   	long_name         atmospheric longwave radiation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FPG                    	long_name         -fraction of potential gpp due to N limitation      units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FPG_P                      	long_name         -fraction of potential gpp due to P limitation      units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    FPI                    	long_name         0fraction of potential immobilization of nitrogen   units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FPI_P                      	long_name         2fraction of potential immobilization of phosphorus     units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FPI_P_vr                      	long_name         2fraction of potential immobilization of phosphorus     units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   FPI_vr                        	long_name         0fraction of potential immobilization of nitrogen   units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   FPSN                   	long_name         photosynthesis     units         umol/m2s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FPSN_WC                    	long_name         Rubisco-limited photosynthesis     units         umol/m2s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FPSN_WJ                    	long_name         RuBP-limited photosynthesis    units         umol/m2s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FPSN_WP                    	long_name         Product-limited photosynthesis     units         umol/m2s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    FROOTC                     	long_name         fine root C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �(   FROOTC_ALLOC                   	long_name         fine root C allocation     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �0   FROOTC_LOSS                    	long_name         fine root C loss   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �8   FROOTN                     	long_name         fine root N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   FROOTP                     	long_name         fine root P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   FROST_TABLE                    	long_name         ,frost table depth (vegetated landunits only)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �P   FSA                    	long_name         absorbed solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �X   FSAT                   	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �`   FSA_R                      	long_name         Rural absorbed solar radiation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �h   FSA_U                      	long_name         Urban absorbed solar radiation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �p   FSDS                   	long_name         $atmospheric incident solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �x   FSDSND                     	long_name         #direct nir incident solar radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSDSNDLN                   	long_name         1direct nir incident solar radiation at local noon      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSDSNI                     	long_name         $diffuse nir incident solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSDSVD                     	long_name         #direct vis incident solar radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSDSVDLN                   	long_name         1direct vis incident solar radiation at local noon      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSDSVI                     	long_name         $diffuse vis incident solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSDSVILN                   	long_name         2diffuse vis incident solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSH                    	long_name         sensible heat      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSH_G                      	long_name         sensible heat from ground      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSH_NODYNLNDUSE                    	long_name         :sensible heat not including correction for land use change     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSH_R                      	long_name         Rural sensible heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSH_U                      	long_name         Urban sensible heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSH_V                      	long_name         sensible heat from veg     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSM                    	long_name         snow melt heat flux    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSM_R                      	long_name         Rural snow melt heat flux      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSM_U                      	long_name         Urban snow melt heat flux      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   FSNO                   	long_name         "fraction of ground covered by snow     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    FSNO_EFF                   	long_name         ,effective fraction of ground covered by snow   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FSR                    	long_name         reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FSRND                      	long_name         $direct nir reflected solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   FSRNDLN                    	long_name         2direct nir reflected solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    FSRNI                      	long_name         %diffuse nir reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �(   FSRVD                      	long_name         $direct vis reflected solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �0   FSRVDLN                    	long_name         2direct vis reflected solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �8   FSRVI                      	long_name         %diffuse vis reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   
F_CO2_SOIL                     	long_name         total soil-atm. CO2 exchange   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   F_CO2_SOIL_vr                         	long_name         0total vertically resolved soil-atm. CO2 exchange   units         gC/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �P   F_DENIT                    	long_name         denitrification flux   units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   
F_DENIT_vr                        	long_name         denitrification flux   units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   F_N2O_DENIT                    	long_name         denitrification N2O flux   units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   	F_N2O_NIT                      	long_name         nitrification N2O flux     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �P   F_NIT                      	long_name         nitrification flux     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �X   F_NIT_vr                      	long_name         nitrification flux     units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �`   GC_HEAT1                   	long_name         #initial gridcell total heat content    units         J/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   GC_ICE1                    	long_name         "initial gridcell total ice content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   GC_LIQ1                    	long_name         "initial gridcell total liq content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   GPP                    	long_name         gross primary production   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   GR                     	long_name         total growth respiration   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   
GROSS_NMIN                     	long_name         gross rate of N mineralization     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    
GROSS_PMIN                     	long_name         gross rate of P mineralization     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   H2OCAN                     	long_name         intercepted water      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   H2OSFC                     	long_name         surface water depth    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   H2OSNO                     	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �    
H2OSNO_TOP                     	long_name         mass of snow in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �(   H2OSOI                        	long_name         0volumetric soil water (vegetated landunits only)   units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �0   H2O_MOSS_WC_TOTAL                      	long_name         total moss water content   units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   HC                     	long_name         heat content of soil/snow/lake     units         MJ/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   HCSOI                      	long_name         soil heat content      units         MJ/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   HEAT_FROM_AC                   	long_name         Lsensible heat flux put into canyon due to heat removed from air conditioning   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   HR                     	long_name         total heterotrophic respiration    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   HR_vr                         	long_name         3total vertically resolved heterotrophic respiration    units         gC/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   HTOP                   	long_name         
canopy top     units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   INT_SNOW                   	long_name         *accumulated swe (vegetated landunits only)     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �P   LABILEP                    	long_name         soil Labile P      units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �X   LABILEP_TO_SECONDP                     	long_name         LABILE P TO SECONDARY MINERAL P    units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �`   
LABILEP_vr                        	long_name         soil labile P (vert. res.)     units         gp/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �h   LAISHA                     	long_name          shaded projected leaf area index   units         none   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LAISUN                     	long_name          sunlit projected leaf area index   units         none   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LAKEICEFRAC                       	long_name         lake layer ice mass fraction   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      P     ��   LAKEICETHICK                   	long_name         @thickness of lake ice (including physical expansion on freezing)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   LAND_UPTAKE                    	long_name         ,NEE minus LAND_USE_FLUX, negative for update   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �H   LAND_USE_FLUX                      	long_name         Atotal C emitted from land cover conversion and wood product pools      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �P   LEAFC                      	long_name         leaf C     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �X   LEAFC_ALLOC                    	long_name         leaf C allocation      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �`   
LEAFC_LOSS                     	long_name         leaf C loss    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �h   LEAFC_TO_LITTER                    	long_name         leaf C litterfall      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �p   LEAFN                      	long_name         leaf N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �x   LEAFP                      	long_name         leaf P     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LEAF_MR                    	long_name         leaf maintenance respiration   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LFC2                   	long_name         3conversion area fraction of BET and BDT that burned    units         per sec    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LITFALL                    	long_name         "litterfall (leaves and fine roots)     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LITHR                      	long_name          litter heterotrophic respiration   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LITR1C                     	long_name         LITR1 C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LITR1C_TO_SOIL1C                   	long_name         !decomp. of litter 1 C to soil 1 C      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   	LITR1C_vr                         	long_name         LITR1 C (vertically resolved)      units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   LITR1N                     	long_name         LITR1 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �0   LITR1N_TNDNCY_VERT_TRANS                      	long_name         -litter 1 N tendency due to vertical transport      units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �8   LITR1N_TO_SOIL1N                   	long_name         !decomp. of litter 1 N to soil 1 N      units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   	LITR1N_vr                         	long_name         LITR1 N (vertically resolved)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   LITR1P                     	long_name         LITR1 P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �0   LITR1P_TNDNCY_VERT_TRANS                      	long_name         -litter 1 P tendency due to vertical transport      units         gP/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �8   LITR1P_TO_SOIL1P                   	long_name         !decomp. of litter 1 P to soil 1 N      units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   	LITR1P_vr                         	long_name         LITR1 P (vertically resolved)      units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   LITR1_HR                   	long_name         Het. Resp. from litter 1   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �0   LITR2C                     	long_name         LITR2 C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �8   LITR2C_TO_SOIL2C                   	long_name         !decomp. of litter 2 C to soil 2 C      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   	LITR2C_vr                         	long_name         LITR2 C (vertically resolved)      units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �H   LITR2N                     	long_name         LITR2 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LITR2N_TNDNCY_VERT_TRANS                      	long_name         -litter 2 N tendency due to vertical transport      units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   LITR2N_TO_SOIL2N                   	long_name         !decomp. of litter 2 N to soil 2 N      units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �@   	LITR2N_vr                         	long_name         LITR2 N (vertically resolved)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �H   LITR2P                     	long_name         LITR2 P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           ��   LITR2P_TNDNCY_VERT_TRANS                      	long_name         -litter 2 P tendency due to vertical transport      units         gP/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     ��   LITR2P_TO_SOIL2P                   	long_name         !decomp. of litter 2 P to soil 2 N      units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            @   	LITR2P_vr                         	long_name         LITR2 P (vertically resolved)      units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x      H   LITR2_HR                   	long_name         Het. Resp. from litter 2   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            �   LITR3C                     	long_name         LITR3 C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            �   LITR3C_TO_SOIL3C                   	long_name         !decomp. of litter 3 C to soil 3 C      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            �   	LITR3C_vr                         	long_name         LITR3 C (vertically resolved)      units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x      �   LITR3N                     	long_name         LITR3 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   LITR3N_TNDNCY_VERT_TRANS                      	long_name         -litter 3 N tendency due to vertical transport      units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     X   LITR3N_TO_SOIL3N                   	long_name         !decomp. of litter 3 N to soil 3 N      units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	LITR3N_vr                         	long_name         LITR3 N (vertically resolved)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   LITR3P                     	long_name         LITR3 P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   LITR3P_TNDNCY_VERT_TRANS                      	long_name         -litter 3 P tendency due to vertical transport      units         gP/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     X   LITR3P_TO_SOIL3P                   	long_name         !decomp. of litter 3 P to soil 3 N      units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	LITR3P_vr                         	long_name         LITR3 P (vertically resolved)      units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   LITR3_HR                   	long_name         Het. Resp. from litter 3   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   LITTERC                    	long_name         litter C   units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           X   
LITTERC_HR                     	long_name         "litter C heterotrophic respiration     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   LITTERC_LOSS                   	long_name         litter C loss      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   
LIVECROOTC                     	long_name         live coarse root C     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   
LIVECROOTN                     	long_name         live coarse root N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   
LIVECROOTP                     	long_name         live coarse root P     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	LIVESTEMC                      	long_name         live stem C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	LIVESTEMN                      	long_name         live stem N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	LIVESTEMP                      	long_name         live stem P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   MR                     	long_name         maintenance respiration    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   M_LITR1C_TO_LEACHING                   	long_name         litter 1 C leaching loss   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   M_LITR2C_TO_LEACHING                   	long_name         litter 2 C leaching loss   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   M_LITR3C_TO_LEACHING                   	long_name         litter 3 C leaching loss   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   M_SOIL1C_TO_LEACHING                   	long_name         soil 1 C leaching loss     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   M_SOIL2C_TO_LEACHING                   	long_name         soil 2 C leaching loss     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   M_SOIL3C_TO_LEACHING                   	long_name         soil 3 C leaching loss     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   M_SOIL4C_TO_LEACHING                   	long_name         soil 4 C leaching loss     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   NBP                    	long_name         Qnet biome production, includes fire, landuse, and harvest flux, positive for sink      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   NDEPLOY                    	long_name         total N deployed in new growth     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   NDEP_TO_SMINN                      	long_name         *atmospheric N deposition to soil mineral N     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   NEE                    	long_name         mnet ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   NEM                    	long_name         DGridcell net adjustment to NEE passed to atm. for methane production   units         gC/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               NEP                    	long_name         Unet ecosystem production, excludes fire, landuse, and harvest flux, positive for sink      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              NET_NMIN                   	long_name         net rate of N mineralization   units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              NET_PMIN                   	long_name         net rate of P mineralization   units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              NFIRE                      	long_name         fire counts valid only in Reg.C    units         counts/km2/sec     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               NFIX_TO_SMINN                      	long_name         1symbiotic/asymbiotic N fixation to soil mineral N      units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   NPP                    	long_name         net primary production     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   OCCLP                      	long_name         soil occluded P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   OCCLP_vr                      	long_name         soil occluded P (vert. res.)   units         gp/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     @   OCDEP                      	long_name         -total OC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   O_SCALAR                      	long_name         8fraction by which decomposition is reduced due to anoxia   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   PARVEGLN                   	long_name         (absorbed par by vegetation at local noon   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   PBOT                   	long_name         atmospheric pressure   units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   PCH4                   	long_name         #atmospheric partial pressure of CH4    units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   PCO2                   	long_name         #atmospheric partial pressure of CO2    units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   PCT_LANDUNIT         
             	long_name         % of each landunit on grid cell    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      H     X   PCT_NAT_PFT                       	long_name         =% of each PFT on the natural vegetation (i.e., soil) landunit      units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �     �   PDEPLOY                    	long_name         total P deployed in new growth     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   PDEP_TO_SMINP                      	long_name         *atmospheric P deposition to soil mineral P     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   PFT_FIRE_CLOSS                     	long_name         Stotal patch-level fire C loss for non-peat fires outside land-type converted region    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   PFT_FIRE_NLOSS                     	long_name         total pft-level fire N loss    units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   PLANT_CALLOC                   	long_name         total allocated C flux     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   PLANT_NDEMAND                      	long_name         &N flux required to support initial GPP     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   PLANT_NDEMAND_COL                      	long_name         &N flux required to support initial GPP     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           X   PLANT_PALLOC                   	long_name         total allocated P flux     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   PLANT_PDEMAND                      	long_name         &P flux required to support initial GPP     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   PLANT_PDEMAND_COL                      	long_name         &P flux required to support initial GPP     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   POTENTIAL_IMMOB                    	long_name         potential N immobilization     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   POTENTIAL_IMMOB_P                      	long_name         potential P immobilization     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   POT_F_DENIT                    	long_name         potential denitrification flux     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	POT_F_NIT                      	long_name         potential nitrification flux   units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   PRIMP                      	long_name         soil primary P     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   PRIMP_TO_LABILEP                   	long_name         PRIMARY MINERAL P TO LABILE P      units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   PRIMP_vr                      	long_name         soil primary P (vert. res.)    units         gp/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   PROD1P_LOSS                    	long_name          loss from 1-yr crop product pool   units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               PSNSHA                     	long_name         shaded leaf photosynthesis     units         umolCO2/m^2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   PSNSHADE_TO_CPOOL                      	long_name         C fixation from shaded canopy      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   PSNSUN                     	long_name         sunlit leaf photosynthesis     units         umolCO2/m^2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   PSNSUN_TO_CPOOL                    	long_name         C fixation from sunlit canopy      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   Q2M                    	long_name         2m specific humidity   units         kg/kg      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   QBOT                   	long_name         atmospheric specific humidity      units         kg/kg      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   QCHARGE                    	long_name         0aquifer recharge rate (vegetated landunits only)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           X   QDRAI                      	long_name         sub-surface drainage   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   QDRAI_PERCH                    	long_name         perched wt drainage    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   QDRAI_XS                   	long_name         saturation excess drainage     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   QDRIP                      	long_name         throughfall    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   QFLOOD                     	long_name         runoff from river flooding     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QFLX_EVAP_TOT                      	long_name         -qflx_evap_soi + qflx_evap_can + qflx_tran_veg      units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QFLX_ICE_DYNBAL                    	long_name         4ice dynamic land cover change conversion runoff flux   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QFLX_LAT_AQU                   	long_name         'Lateral flow between hummock and hollow    units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QFLX_LIQ_DYNBAL                    	long_name         4liq dynamic land cover change conversion runoff flux   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QFLX_SURF_INPUT                    	long_name         Runoff from hummock to hollow      units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QH2OSFC                    	long_name         surface water runoff   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QINFL                      	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QINTR                      	long_name         interception   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QIRRIG                     	long_name         water added through irrigation     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QOVER                      	long_name         surface runoff     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	QOVER_LAG                      	long_name         +time-lagged surface runoff for soil columns    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QRGWL                      	long_name         9surface runoff at glaciers (liquid only), wetlands, lakes      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QRUNOFF                    	long_name         0total liquid runoff (does not include QSNWCPICE)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   QRUNOFF_NODYNLNDUSE                    	long_name         ]total liquid runoff (does not include QSNWCPICE) not including correction for land use change      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	QRUNOFF_R                      	long_name         Rural total runoff     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	QRUNOFF_U                      	long_name         Urban total runoff     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               QSNOMELT                   	long_name         	snow melt      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              	QSNWCPICE                      	long_name         #excess snowfall due to snow capping    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              QSNWCPICE_NODYNLNDUSE                      	long_name         Pexcess snowfall due to snow capping not including correction for land use change   units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              QSOIL                      	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               QVEGE                      	long_name         canopy evaporation     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   QVEGT                      	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   RAIN                   	long_name         atmospheric rain   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   RETRANSN                   	long_name         plant pool of retranslocated N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   RETRANSN_TO_NPOOL                      	long_name         deployment of retranslocated N     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   RETRANSP                   	long_name         plant pool of retranslocated P     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   RETRANSP_TO_PPOOL                      	long_name         deployment of retranslocated P     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           X   RH2M                   	long_name         2m relative humidity   units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   RH2M_R                     	long_name         Rural 2m specific humidity     units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   RH2M_U                     	long_name         Urban 2m relative humidity     units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   RR                     	long_name         /root respiration (fine root MR + total root GR)    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   SABG                   	long_name         solar rad absorbed by ground   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SABG_PEN                   	long_name         2Rural solar rad penetrating top soil or snow layer     units         watt/m^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SABV                   	long_name         solar rad absorbed by veg      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SCALARAVG_vr                      	long_name         average of decomposition scalar    units         fraction   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SECONDP                    	long_name         soil secondary P   units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	   SECONDP_TO_LABILEP                     	long_name         SECONDARY MINERAL P TO LABILE P    units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	   SECONDP_TO_OCCLP                   	long_name         !SECONDARY MINERAL P TO OCCLUDED P      units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	    
SECONDP_vr                        	long_name         soil secondary P (vert. res.)      units         gp/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     	(   	SEEDC_GRC                      	long_name         /pool for seeding new PFTs via dynamic landcover    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN                      	long_name         soil mineral N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_NPOOL                     	long_name         #deployment of soil mineral N uptake    units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_PLANT                     	long_name         plant uptake of soil mineral N     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_SOIL1N_L1                     	long_name         +mineral N flux for decomp. of LITR1to SOIL1    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_SOIL2N_L2                     	long_name         +mineral N flux for decomp. of LITR2to SOIL2    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_SOIL2N_S1                     	long_name         +mineral N flux for decomp. of SOIL1to SOIL2    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_SOIL3N_L3                     	long_name         +mineral N flux for decomp. of LITR3to SOIL3    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_SOIL3N_S2                     	long_name         +mineral N flux for decomp. of SOIL2to SOIL3    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINN_TO_SOIL4N_S3                     	long_name         +mineral N flux for decomp. of SOIL3to SOIL4    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINP                      	long_name         soil mineral P     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINP_LEACHED                      	long_name         $soil mineral P pool loss to leaching   units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           	�   SMINP_TO_PLANT                     	long_name         plant uptake of soil mineral P     units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
    SMINP_TO_PPOOL                     	long_name         #deployment of soil mineral P uptake    units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
   SMINP_TO_SOIL1P_L1                     	long_name         +mineral P flux for decomp. of LITR1to SOIL1    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
   SMINP_TO_SOIL2P_L2                     	long_name         +mineral P flux for decomp. of LITR2to SOIL2    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
   SMINP_TO_SOIL2P_S1                     	long_name         +mineral P flux for decomp. of SOIL1to SOIL2    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
    SMINP_TO_SOIL3P_L3                     	long_name         +mineral P flux for decomp. of LITR3to SOIL3    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
(   SMINP_TO_SOIL3P_S2                     	long_name         +mineral P flux for decomp. of SOIL2to SOIL3    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
0   SMINP_TO_SOIL4P_S3                     	long_name         +mineral P flux for decomp. of SOIL3to SOIL4    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
8   SMINP_vr                      	long_name         soil mineral P (vert. res.)    units         gp/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     
@   SMIN_NH4                   	long_name         soil mineral NH4   units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           
�   SMIN_NH4_vr                       	long_name         soil mineral NH4 (vert. res.)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     
�   SMIN_NO3                   	long_name         soil mineral NO3   units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   SMIN_NO3_LEACHED                   	long_name         soil NO3 pool loss to leaching     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   SMIN_NO3_RUNOFF                    	long_name         soil NO3 pool loss to runoff   units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   SMIN_NO3_vr                       	long_name         soil mineral NO3 (vert. res.)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     P   SNOBCMCL                   	long_name         mass of BC in snow column      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SNOBCMSL                   	long_name         mass of BC in top snow layer   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SNODSTMCL                      	long_name         mass of dust in snow column    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SNODSTMSL                      	long_name         mass of dust in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SNOINTABS                      	long_name         7Percent of incoming solar absorbed by lower snow layers    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SNOOCMCL                   	long_name         mass of OC in snow column      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SNOOCMSL                   	long_name         mass of OC in top snow layer   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SNOW                   	long_name         atmospheric snow   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               SNOWDP                     	long_name         gridcell mean snow height      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              SNOWICE                    	long_name         snow ice   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              SNOWLIQ                    	long_name         snow liquid water      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              
SNOW_DEPTH                     	long_name          snow height of snow covered area   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               
SNOW_SINKS                     	long_name         snow sinks (liquid water)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   SNOW_SOURCES                   	long_name         snow sources (liquid water)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   SOIL1C                     	long_name         SOIL1 C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   SOIL1C_TO_SOIL2C                   	long_name         decomp. of soil 1 C to soil 2 C    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   	SOIL1C_vr                         	long_name         SOIL1 C (vertically resolved)      units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     H   SOIL1N                     	long_name         SOIL1 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOIL1N_TNDNCY_VERT_TRANS                      	long_name         +soil 1 N tendency due to vertical transport    units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL1N_TO_SOIL2N                   	long_name         decomp. of soil 1 N to soil 2 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   	SOIL1N_vr                         	long_name         SOIL1 N (vertically resolved)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     H   SOIL1P                     	long_name         SOIL1 P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOIL1P_TNDNCY_VERT_TRANS                      	long_name         +soil 1 P tendency due to vertical transport    units         gP/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL1P_TO_SOIL2P                   	long_name         decomp. of soil 1 P to soil 2 N    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   	SOIL1P_vr                         	long_name         SOIL1 P (vertically resolved)      units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     H   SOIL1_HR                   	long_name         Het. Resp. from soil 1     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOIL2C                     	long_name         SOIL2 C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOIL2C_TO_SOIL3C                   	long_name         decomp. of soil 2 C to soil 3 C    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SOIL2C_vr                         	long_name         SOIL2 C (vertically resolved)      units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL2N                     	long_name         SOIL2 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   SOIL2N_TNDNCY_VERT_TRANS                      	long_name         +soil 2 N tendency due to vertical transport    units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     X   SOIL2N_TO_SOIL3N                   	long_name         decomp. of soil 2 N to soil 3 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SOIL2N_vr                         	long_name         SOIL2 N (vertically resolved)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL2P                     	long_name         SOIL2 P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   SOIL2P_TNDNCY_VERT_TRANS                      	long_name         +soil 2 P tendency due to vertical transport    units         gP/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     X   SOIL2P_TO_SOIL3P                   	long_name         decomp. of soil 2 P to soil 3 N    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SOIL2P_vr                         	long_name         SOIL2 P (vertically resolved)      units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL2_HR                   	long_name         Het. Resp. from soil 2     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   SOIL3C                     	long_name         SOIL3 C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           X   SOIL3C_TO_SOIL4C                   	long_name         decomp. of soil 3 C to soil 4 C    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   	SOIL3C_vr                         	long_name         SOIL3 C (vertically resolved)      units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     h   SOIL3N                     	long_name         SOIL3 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOIL3N_TNDNCY_VERT_TRANS                      	long_name         +soil 3 N tendency due to vertical transport    units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL3N_TO_SOIL4N                   	long_name         decomp. of soil 3 N to soil 4 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   	SOIL3N_vr                         	long_name         SOIL3 N (vertically resolved)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     h   SOIL3P                     	long_name         SOIL3 P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOIL3P_TNDNCY_VERT_TRANS                      	long_name         +soil 3 P tendency due to vertical transport    units         gP/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL3P_TO_SOIL4P                   	long_name         decomp. of soil 3 P to soil 4 N    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   	SOIL3P_vr                         	long_name         SOIL3 P (vertically resolved)      units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     h   SOIL3_HR                   	long_name         Het. Resp. from soil 3     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOIL4C                     	long_name         SOIL4 C    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SOIL4C_vr                         	long_name         SOIL4 C (vertically resolved)      units         gC/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL4N                     	long_name         SOIL4 N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   SOIL4N_TNDNCY_VERT_TRANS                      	long_name         +soil 4 N tendency due to vertical transport    units         gN/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     p   SOIL4N_TO_SMINN                    	long_name         #mineral N flux for decomp. of SOIL4    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SOIL4N_vr                         	long_name         SOIL4 N (vertically resolved)      units         gN/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL4P                     	long_name         SOIL4 P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   SOIL4P_TNDNCY_VERT_TRANS                      	long_name         +soil 4 P tendency due to vertical transport    units         gP/m^3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     p   SOIL4P_TO_SMINP                    	long_name         #mineral P flux for decomp. of SOIL4    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SOIL4P_vr                         	long_name         SOIL4 P (vertically resolved)      units         gP/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOIL4_HR                   	long_name         Het. Resp. from soil 4     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   SOILC                      	long_name         soil C     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   SOILC_HR                   	long_name          soil C heterotrophic respiration   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   
SOILC_LOSS                     	long_name         soil C loss    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOILICE                       	long_name         #soil ice (vegetated landunits only)    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   SOILLIQ                       	long_name         ,soil liquid water (vegetated landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x         SOILPSI                       	long_name         'soil water potential in each soil layer    units         MPa    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     x   SOILWATER_10CM                     	long_name         @soil liquid water + ice in top 10cm of soil (veg landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SOLUTIONP                      	long_name         soil solution P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SOLUTIONP_vr                      	long_name         soil solution P (vert. res.)   units         gp/m^3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x         SOMHR                      	long_name         -soil organic matter heterotrophic respiration      units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   SOM_C_LEACHED                      	long_name         .total flux of C from SOM pools due to leaching     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SR                     	long_name         'total soil respiration (HR + root resp)    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   STORVEGC                   	long_name         )stored vegetation carbon, excluding cpool      units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   STORVEGN                   	long_name         stored vegetation nitrogen     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   STORVEGP                   	long_name         stored vegetation phosphorus   units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SUPPLEMENT_TO_SMINN                    	long_name         supplemental N supply      units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SUPPLEMENT_TO_SMINP                    	long_name         supplemental P supply      units         gP/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   	SoilAlpha                      	long_name         factor limiting ground evap    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   SoilAlpha_U                    	long_name         !urban factor limiting ground evap      units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TAUX                   	long_name         zonal surface stress   units         kg/m/s^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TAUY                   	long_name         meridional surface stress      units         kg/m/s^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TBOT                   	long_name         atmospheric air temperature    units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TBUILD                     	long_name         #internal urban building temperature    units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TG                     	long_name         ground temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TG_R                   	long_name         Rural ground temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TG_U                   	long_name         Urban ground temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TH2OSFC                    	long_name         surface water temperature      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               THBOT                      	long_name         %atmospheric air potential temperature      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              TKE1                   	long_name         (top lake level eddy thermal conductivity   units         W/(mK)     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              TLAI                   	long_name         total projected leaf area index    units         none   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              TLAKE                         	long_name         lake temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      P         TOTCOLC                    	long_name         >total column carbon, incl veg and cpool but excl product pools     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   	TOTCOLCH4                      	long_name         9total belowground CH4, (0 for non-lake special landunits)      units         gC/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   TOTCOLN                    	long_name         +total column-level N but excl product pools    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTCOLP                    	long_name         +total column-level P but excl product pools    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   
TOTECOSYSC                     	long_name         Ftotal ecosystem carbon, incl veg but excl cpool but excl product pools     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   
TOTECOSYSN                     	long_name         (total ecosystem N but excl product pools   units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   
TOTECOSYSP                     	long_name         (total ecosystem P but excl product pools   units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTLITC                    	long_name         total litter carbon    units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   
TOTLITC_1m                     	long_name         $total litter carbon to 1 meter depth   units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTLITN                    	long_name         total litter N     units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTLITP                    	long_name         total litter P     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   
TOTLITP_1m                     	long_name         total litter P to 1 meter      units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTPFTC                    	long_name         )total patch-level carbon, including cpool      units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTPFTN                    	long_name         total PFT-level nitrogen   units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTPFTP                    	long_name         total PFT-level phosphorus     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTSOMC                    	long_name          total soil organic matter carbon   units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   
TOTSOMC_1m                     	long_name         1total soil organic matter carbon to 1 meter depth      units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTSOMN                    	long_name         total soil organic matter N    units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TOTSOMP                    	long_name         total soil organic matter P    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               
TOTSOMP_1m                     	long_name         &total soil organic matter P to 1 meter     units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              TOTVEGC                    	long_name         (total vegetation carbon, excluding cpool   units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              TOTVEGC_ABG                    	long_name         4total aboveground vegetation carbon, excluding cpool   units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              TOTVEGN                    	long_name         total vegetation nitrogen      units         gN/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               TOTVEGP                    	long_name         total vegetation phosphorus    units         gP/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   TREFMNAV                   	long_name         (daily minimum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   
TREFMNAV_R                     	long_name         .Rural daily minimum of average 2-m temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   
TREFMNAV_U                     	long_name         .Urban daily minimum of average 2-m temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   TREFMXAV                   	long_name         (daily maximum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   
TREFMXAV_R                     	long_name         .Rural daily maximum of average 2-m temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   
TREFMXAV_U                     	long_name         .Urban daily maximum of average 2-m temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           X   TSA                    	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   TSAI                   	long_name         total projected stem area index    units         none   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   TSA_R                      	long_name         Rural 2m air temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   TSA_U                      	long_name         Urban 2m air temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   TSOI                      	long_name         +soil temperature (vegetated landunits only)    units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   	TSOI_10CM                      	long_name         $soil temperature in top 10cm of soil   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TSOI_ICE                      	long_name         %soil temperature (ice landunits only)      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x         TV                     	long_name         vegetation temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   TWS                    	long_name         total water storage    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   TWS_MONTH_BEGIN                    	long_name         /total water storage at the beginning of a month    units         mm     cell_methods      time: instantaneous    
_FillValue        {@��   missing_value         {@��           �   TWS_MONTH_END                      	long_name         )total water storage at the end of a month      units         mm     cell_methods      time: instantaneous    
_FillValue        {@��   missing_value         {@��           �   T_SCALAR                      	long_name         'temperature inhibition of decomposition    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   U10                    	long_name         	10-m wind      units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              URBAN_AC                   	long_name         urban air conditioning flux    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              
URBAN_HEAT                     	long_name         urban heating flux     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               VOLR                   	long_name         !river channel total water storage      units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   VOLRMCH                    	long_name         (river channel main channel water storage   units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   WA                     	long_name         :water in the unconfined aquifer (vegetated landunits only)     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   	WASTEHEAT                      	long_name         Csensible heat flux from heating/cooling sources of urban waste heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   WF                     	long_name         )soil water as frac. of whc for top 0.05 m      units         
proportion     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   WIND                   	long_name         #atmospheric wind velocity magnitude    units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           P   WOODC                      	long_name         wood C     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           X   WOODC_ALLOC                    	long_name         wood C eallocation     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           `   
WOODC_LOSS                     	long_name         wood C loss    units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           h   WOOD_HARVESTC                      	long_name         &wood harvest carbon (to product pools)     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           p   WOOD_HARVESTN                      	long_name         !wood harvest N (to product pools)      units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           x   WTGQ                   	long_name         surface tracer conductance     units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �   W_SCALAR                      	long_name         .Moisture (dryness) inhibition of decomposition     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     �   XR                     	long_name         total excess respiration   units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               XSMRPOOL                   	long_name         temporary photosynthate C pool     units         gC/m^2     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              ZBOT                   	long_name         atmospheric reference height   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              ZWT                    	long_name         ,water table depth (vegetated landunits only)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��              ZWT_CH4_UNSAT                      	long_name         Fdepth of water table for methane production used in non-inundated area     units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��               	ZWT_PERCH                      	long_name         4perched water table depth (vegetated landunits only)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           (   	cn_scalar                      	long_name         N limitation factor    units             cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           0   	cp_scalar                      	long_name         P limitation factor    units             cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           8   leaf_npimbalance                   	long_name         9leaf np imbalance partial C partial P/partial C partial N      units         gN/gP      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           @   nlim_m                     	long_name         runmean N limitation factor    units             cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           H   o2_decomp_depth_unsat                         	long_name         o2_decomp_depth_unsat      units         mol/m3/2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      x     P   plim_m                     	long_name         runmean P limitation factor    units             cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��           �;�r<���=�=�o�>YI:>�l�?�~?��?�#'@7U�@��,@��rAN�A���B��=L��?��@ff@�33A��AI��A���A���B	L�B3�;�r<���=�=�o�>YI:>�l�?�~?��?�#'@7U�@��,@��rAN�A���B��º�Bº�BB>@�B>@�EA4�EA4�        ?�  ?�              G�  
%      ��     )�@��     @�7�    11/27/23        16:13:01        4�G4ݒ�16�\1S{��f��h�6� �6���5��6�BB��B��B��B��B��B��7��7�s�7���7�F61�h12r�                *�ۥ*�ۥ6��6�#�        ?Eo*?C�G        4)��53Lp3�43�I�2�3,S�0�i0<}�/�JT1�'
2#�1$�C    0/0�        =��<=.�$>8!=�"B>V��=��>5\=|��>6i=��=���=���<�}�=Nk�<	t<���;_�;�?�:��:�"                                        8C �;�;8C�q<)8�9w<̡�9ɩ`=8?�<�^�=�Y�=IA�=�e=)x�=tܚ<�D�<ӯ�;�O;�`o:�_:��L                                        >�j>�Wy=�y=��L=�M�=��W=��L=S!x=�g<�4;]5O6S��96� 4|7mo�2h�	.y��+/~�+�+�                                        @��N@U�@�}@)#z@�n�?�3�@7�P?S7�?���>��=պ�7�]�9��u1��1���+tS+$�+!�+�+�                                        DC���C1��C.�T        60�6)~�5��Z6 �/5�,5"tqD�D���D���D_ǆD�D�XC�b�C��|C}BCT�B�YC��Bi��BQ�AďAϢ>��>��                                                >�`>�Q�1x��1� �0� 90���@~�I@;��@�v?���?�NU?�=?	��?�>�*�>ާ>ie(>�7�=�A�=�r�<�`D<�s:E�:Y�                                                =t��=ő0"��00a�/M޴/^�3?'p>���>�&�>�l�>FO�>1�8=��/=�1�=TQ|=�q�=Dl=ON�<�^�<��";QJ;D�58���8�L�4��w4��U                                        B�{�B��>
>'�<��<�>�C�|BC�6|?��?"�=��M=��2��2��(�/�<IjD,ID7�nA��A��?6��??S/X�N/X�N        2d��u*{@��{@��                                                        {@��{@��B!{oB]}B!{oB]}{@��{@��@�F@�b8)�8�$e�����D!��"���� �{�(�4.��.�����]����?,��?*P�-ÅU-L*9@I��@Is.��/��:3e�4�a��=��@Y>�ci?��A��AѸ�A#��A�?�p?��=u�>=���?�p?��{@��{@��9K&>a�49K&>a�49K&>a��A�2�A�OFA�2�A�OF{@��{@��C��C���C��C���{@��{@��C�	�C�	�?xFT?|OY?�  ?�  ?8+�?I�x?N��?a]�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�                                          ?�?3zI?6��?P՛?W�?o�?l�R?| �?x��?|�[?z$9?}u?z|5?}~?{�?}�q?{Q+?}�"?{\+?}��                                        @-�@+O?��?�%R?�*�?��=��u=��iC-wcC3�6�&�6��6�)h6��@�*3@��t>1�*>7�?�c?ĭ�B��B���=�}>�D�B��B���{@��{@��Bݗ-Bݗ-B+�B+�C[hC[hAHX�AHX�B��B��B�PGB�PGA���A���Bo�-Bo�-B�B	Պ@L�[@_��B�B	ՊB�B	Պ{@��{@��A��5A��Y@��@I.+@��@I.+{@��{@��>�9>�q�>�,�>��7A��A�AfA��A�4A�OTA���@HE�@J��@8��@HjA.�iA<ֆ@��@(�@                                                                                                                                2��2��(2���3���31�,3���3��)3��R4XI3���3''2��12P31Z6/��/���*�4|(WB(I��(���&�'h                                        0�0g�-{�#-0�=2���2��|2���3�Ψ35��3�93䞋3��4�3��E3�F2��42i�1Z��0'7�0���*�J�(o~(I��(��R&��<'}                                        N�nfN�)B��uB�u�E�7�E�8�c8��6���6��50_k53�;2���2���=0_=,/�:�%@��UA�&4A��?q��?��>���?��>�0?&�>��p?9 ?��?Q��?Cj�?[��?TKo?X�$?Qr�?R��?R,t?SM?t<?{y>��L>���                                        @z@��IFc�FdآFc��Fd�,        7cnt7ig�99#�8�m8ᮽ8���8�np8}v8.7o��7_�o6��6p�z5=� 5c{$3�oV4
1}1�2-L/$=�/��                                        A%q�A*��C�_�C��>�9{>�Ss0b��0Z]�@l��@��0@IG�@I'f@�<?��e?��x?���>�u�>��>(=̡<Ҡ�<��y;�G�;��9���: 77�D7Bh                                        ?�)�?�{h>��>��q{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@�γo\�ΐ        C�CQ�6��56�/6���6�"6���6���@��@��V>���>ܨ�6���6�W�        7^<�7hM�6�$�6�(�A$��@���6
�6r>CR�3C=rB�.B�diB�?A���A�#@� �@@��@���?��'@P?��?�=�s�=䲷<�:1<���:L/):L)�                                        >z�	>Pl���A����3�z2��2n5V2l��1'*�0�n/.ǲ.{ײ-!>�-�$I,�4�,߹7+(�+:Kp)X#�)a�'x�'�Kz                                        3��3� *@�ҙ@y�#@0Wu@
��?q�:?x->Wb�>W��=���> ��=!�=_[�<M�<_�:;%F�;,GB9�u�9�(?7�l7��                                        <�mh<f����*�Щ�180�1	z0��0��n/6�.�]�-8s,�w+C��+�f�*]�*��l(��)R�'�'��%'<�%-�                                        1���1źG>�|�>�y�>L~)>!�g=��==&��<U<[��;�_;��:�5�;*�:V�:��8�s�8��7��z7�l5P5O�B                                        5��q5�x�B�	yB�U�6T�'6^��D��]Dt�jD=  D0U�C��C�+�C��C8�iB�ߩB��AB&v�B�;eA��A��v@z�@��>���>�e;��];�o�                                        ?�a@���cƴWC3i<2��3PD3��2A��2=��0�Q�0��S/OU�/^�~.�&/!{-�,b-���+�R�+���)��)�}                                        3���4�)A��cA���Ab�AV+Z@��t@��n@B�_@l�j?�8@]D?Th�?��x>�I>>ʹ�=��=�uH<��<�59u�9��                                        =��e>	?{�ۧ[��{1�{1�iU1���1gTN0i��0c�.�.��-J8-i(,[�,-o�+t��+��J)��)�X�'���'��	                                        2Ur2!�?�%�?�'�?��"?��>��D? �>;�>e&�=�#C>��=7�=�;�<���<��;�z;�Iq9�[39�)�6ϵ�6�L�                                        6��6�-�C8�C/6V�6vD���D�!Dq�jDR
�D5�D�_C{�`C��lCI{Ck�B���C��B/� BQ�A�tA�>�Z? z�;�+t;��\                                        @9RP@\Ah�a/Q���2;A�2 �Q37ښ2�q2wǔ2@��1��0SB\/`�.�B.��H/�Hk.5H�.c*�,�s�,�i�)��)��?                                        3{�$3��A��wA��6A�(�Ar�YA!{(A"�@�3@�OY@>Y�@�M�?��@T�L?jl,?��>L!>W�<:J�<<�9 ys8���                                        >>z�>Y�u����u 0�/0���1�tO19�0���0��3/ 3.�k#-��-SjB,�e�-�N�,Z,<Y�*n�^*��&'�"'�y                                        1��1��@?��?��?��?�!?6�?6<>��m>Ɣ*>'�>��=��>0�f=;�=`u<'[<"�:>�:�6���6���                                        5x�5�r�C~�C���6�$�6�(�7dau7k��@g$@�6=*�u=:M:6#�:F��A�SA��>&8�>5Z�;1M�;Aq�7@e7F4�                                                        3o\4ΐ4�Ϲ4��V1zY;1zY;�o\�ΐ��"���~3��4&4�w�4��2�Ag2�0�e0�2�2�V�7es7k؀::Ͳ:-B;��;��;�R�;��;eX�;U��:��u:��:Cʈ:��9Ay9FK7�r~8��6|�r6��K4d��4��1'�#1�ܒ                                        ,���,���?�?y�?�  ?c��?�  ?K�?��?&�@?@�'>�>���>Q�E>b�>N�	>R�a>N��>MI�>L��>L��>L��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��CXGCOgG�a\G�a\{@��{@��A۲A۲B�  B�                                                                                  B  B  Ap  Ap                                                          A�  A�  A�  A�                                  2��32՟�-�Ȩ-�Ȩ32�x0;�/��7��	7�JT4��4�]�4��4�]�2�[Z2�̫2��c2��y2��c2��y5aE*5,;�1���1���2�A2�J�2͇Q2���<��<��-F�d-F�d?F?F>���>���>;�>;�=6%=6%<-Ȯ<-Ȯ:��:��8���8���5$�_5$�_/���/���'=�D'=�D                                                ?C�?�7��M7��?��?�� 7�-y7�G;�s�;�N�;�.;�.8	r:x�6{К6v�|                7�7�7�B�        7�Z�7��,        ����6J��            6��1X7^aV7��e9�3�6���6���        6D��    6D��            6�U�7�[6�U�7�[6�U�7�[{@��{@��6�+�7�&    H��    H��6��6���5�>h5���7597/ۀ7��;7��;@��%@�)3�-�3�}�>���>��g1��c1��;B�Z�B�1B�Z�B�1{@��{@��7�x7"A�\?A��z?�k?��eB���B��>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�?ј�?��&/�v/�V+1Ѫ+&*IAj�hAu�~AD �A?6@�#m@�@�9X@��?��?�H�>�R>��)=��f=��?<H@�<s�:S��:��7Kp�7�V�                                                >��%?D��4��N4�!�4��N4�!�3�J�3��%4&N4%(X�N4r�K�c44�x4<�?���"�-Q����@�?��:+t��+a�2�+Z2���2�+Z2���/����y�K���D[=����X�/�/���β	���㱹m�A���A�Avb?Aql�A ��A��@���@�w�@w�?��?��>��%=�=�=�	<���<��V:�u�:��L7���7�k                                        >�g$?B�{>~q?Y��>u��?�z^?5�?�	?��@7r?R�?�A�>���?%
�>D�>��I=?�=�i;�}Q;�M;8�z98��                                        <� <�.D�-�� /8�    9'm7:�2�9��
:���;�-�;���=Hb5<���=v�=��=x)�;��:���;t|0x6O.��.1%_.��,�,�C)                                        4}��4�ѡ3W�3<u�:2U�:N�L8�o?9$�=/�s=)��6��6��4��4�P�6�zS6�zS=�A@=��A�}%A�1E>yU�>�\\=�Xh=�l�6��6�[6խ|6���A�ˢA��}5��%5�|iC(�C>�BպKB��~B��B��B�BU�&A��`B�mA9�JA�_#@�D�@��?��F?���>��>�U                                                ?�d�@�S��>�@A@1���KH\2���1s�1䃦2{@0�E�0��/~�/��Y.�6�/AQ�-��J-���,�	,�)ߓ�)�d                                        48�45�tA-��A4S�A|2A#�T@�mW@��@K�n@�X�?�3*@3|�?w��?�)�>Ͱ�>�+�=�G=�-*</��<1q                                                =S��=��m�AT.��g.����շ0	��/��/s��/���.<d.:H�-| -�_,$~�,�4�+(KY+EL)�v�)�@�'k�9'm�:                                        1�b�1�k>��>�YG>���>���>F�>��u=�v�>�&=qG�=�s�=�=\�<[gn<x��;F�2;RRU9�n�9���6��[6�-!                                        5�V5ҊC�CN�h6(��6,�D^��DV�qDCbDD>X_DtBD�IC���C���C7��C��}B��&CZرBq�<B�
�AD��ARvq>��>�Z;5�; 1                                        AS�A��E����� ��,���3�!3��#3p:�37�*2��1@�.0H`/p�0A�0�/��/�-�-�g-�*��*��s                                        4�@X4�CB���B�LB�A�B}�)B0�B?�bA���BآAt�%AΔ�A(p�A���@�2�@�n?�?�N�=&h�=(��90G�9.՗                                        >��d?�c;[�I���cc*��M�1J��1��1 S0��/���.�X�-ջ�,�VZ-��w.w
�-	�-00�+b��+sQ�(`J6(cex                                        2^"�2cz�@i/@��@
��@[I?�a?̄�?LI[?���?�i?\ZK>��{?��>+��>Q u=��=��:���:��7+�x7*^Z                                        6��6&�C|�C�7�6b�6(�DP,�D%�DG�;D0�D*UzD;^C���C��C��	C�k'CD��Cs`kBߓqB�mB8��BQy�A	�AD=�>�b?V�q                                        Aɟ Aӌ����}�ie�������!2c4���4k�U4�3�܇2���2���1��2y��1#�S1y�107�=0w.�n�.�o?                                        5 W�5�kB��%B��B���B~�gB�DaBa��B<ĉB3'JA��tB	�SA��
A³�A2�ZAI�W@��@��h?\a�?��N=2��=��'                                        ?	�?c������A貕_��Ls�����/J�-1�>�1���1%$�0��[/��/�Ph///��X.Q�Y.��F-kr-�a+�Y�, F�                                        2$GZ2(�t?�+�?�	�?̵|?��?�l?��?q�g?eQ?<s?0�~>ɱ�>�7�>d�>��=�'�=ր�<�%<��O:d�:�*                                        60x�65SC˼PC���Ce=�CSl�Cc��CR�MC_"cCO>|CR�7CG�C>�vC:�#C'��C)�jC)wC��B�VB��QB� BԐB�ǲB�jy                                        B"�B.]F�'����V��t��L:��`3�`�Y��)�0n3�33F%U3`D	3v�2�j�2�%1�x21]u1GD�1HO�1Sϣ                                        4`�4d�iA�dgA�#�A�b*A�j�A���A�˗A��,A��A��+A���A�O(A���AhBWAm�#ADD�AK��A��A*�@��A��                                        ?P�?_/��V��Kic�5jѲ�� A���B��S�-(<1r�0}�o0���0.��/ӻ�0��/h�/@.;X�.N.�2�.��                                        1��01�K%>��>�~>�si>ג�>�}S>�7�>��>�a\>�NL>�x�>��q>��>��X>�i>{9l>�]Z>I�j>Y�>.�>,u@                                        6�G6�aDQ8�Dh��6��#6���6��#6���?Ճn@:�[@8�]@���@�I%@��A
�Ad��Ay�6A��At�A���@��J@��                                                                @v�|@���@�b�AS��Aw�aA��Br�B7��B��\B��C[�C�C�n�C�C�c�C�W�D��D�6D7D7|                                        �W�$�Q%~�Y�s�\���IN�U��(�<�9X�����	�+Ed�M R�Ѳ�>�����F[��9Y��
���v���p  �p  �p  �p  �p  �p  �p  �p  �p  �p  B'�TB���:��:�];��6;��T;x�;x�-;"p�;z�:���:�TG:9�@"9"o�8�x8jy8]�6���6��
4��5�2!�2n}n                                        6��#6���        7���7ż�D�G�Dk�B�&B'�@Nw=@T��        0m�0�J�?�,?��{@��{@�μms-�        C�uC�u{@��{@��C��hC���C��%C��{@��{@��C�,^C�W$C�uC�u{@��{@��@�8@C4{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@���  �  ?�?q�{�  �  �  �  EA�E?�C��C� A�]A
^C~�C���Cx&5C�s @��	@��">�r@>��O>�b�>�n!D�Z�D�vB@eOBM��@|4�@�.\DQ8�Dh��DUhD*$B�v�B��?��@�
?�˔?�T�D�I}D�g�D)�D4� B@eOBM��@|4�@�.\C���C���C���C���{@��{@��C�0�C�)jC�0�C�)j{@��{@��C��uC��?X	�?[�sC��uC��{@��{@��C�'�C�6IC�/C�02C�
cC� �C��C�	dC�ەC���C�ŊC���C���C��TC���C��WC���C���C���C���C��vC��GC��TC�¬C��&C��JC���C��*C���C��C�YC�&�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C�/E�sZE��{@��{@��E��]E�f�?C?��?|?o�?̩?�Q? ��?��>��>�خ>�~�>��>�{�>��L>�d�>�g>�m>��>��>��3{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��?| N?~�                                E�@ E�@         >���?,�?b�?b�C���C̎p68�6*�,6L�6*w�                ;��:�>�?-e?ޝ?X�?$a?M#?,�}? �?5��?6SK?F�E?Zj�?]�?t�?s&�?}�?~5�?�w?�?�\?�w{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��6hF�6{��]׎�k��A�  A�  >v+_=��>�� >
`?>?�q�{@��{@��{@��{@��        {@��{@��7m�c7�v7D7H֞7��6�}�6��6�>�5���5��t4�Z4n�3���3�0�2��c2�q�1U�_1X4/��*/�~x                                        {@��{@��