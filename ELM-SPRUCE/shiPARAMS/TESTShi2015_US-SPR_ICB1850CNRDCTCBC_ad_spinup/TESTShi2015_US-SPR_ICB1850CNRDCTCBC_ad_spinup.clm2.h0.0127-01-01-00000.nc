CDF      
      lndgrid       gridcell      landunit      column         pft    @   levgrnd       levurb        levlak     
   numrad        levsno        ltype      	   natpft        string_length         levdcmp       levtrc     
   hist_interval         time             title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 11/27/23 16:21:18   source        Community Land Model CLM4.0    hostname      ubuntu     username      mwjones    version       5823c39    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
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
_FillValue        {@��   missing_value         {@��           �;�r<���=�=�o�>YI:>�l�?�~?��?�#'@7U�@��,@��rAN�A���B��=L��?��@ff@�33A��AI��A���A���B	L�B3�;�r<���=�=�o�>YI:>�l�?�~?��?�#'@7U�@��,@��rAN�A���B��º�Bº�BB>@�B>@�EA4�EA4�        ?�  ?�              G3�  aU      ��     א@�5`    @�t�    11/27/23        16:21:18        4�2�4��914�w1W�F.��o-��6�P�6��5��6
IB��B��B��B��B��B��7�ƥ7�17��7�*�1��*2��                *�ۥ*�ۥ6�086��F        ?EQc?C��        4'NN5*f�3���3�:72�3'a�0�E�0I��/���1��1�o1#�    0ԥ        =��=&��>3�a=��)>W�Z=~�0>3��=oy�>�|=�II=�"�=�<<�t�=N�@<)ݮ<���;�g[;���:�,�;Mh                                        8C��;�6�8H<�<Zo8���<���9�0C=+n<]=��{=L	=ӊ�=-��=xZL<��.<��;Β/<X�:�E9; �                                        >�W�>��p=��==�D=��=�а=bȦ=��<�;=t�710x9l4���6�6n2[�-� [+.�+�+�                                        @�-�@X�T@���@,�$@�|�?���@;��?X�?�e�>('�=�=�8� 9���1��1]2V+t��+##G+!�+�+�                                        DF�C��C2IIC2�]        6�x6,�5�X6�X54�5%b�D�N�D�جD�/�Dc��D�fD�.C���C�eCJCS�QB��WC�iBrnEB_��A?F�A7�g>�-t>��f                                                >�0`>��x1t�o1�q�0��F0���@|S)@A{@:?�
�?��t?�V	?	?.�>��#>�{�>h�>�Y�=��-=���<��<�:�:Z�5:c�!                                                =vE=u��0 f�03�/J��/b�,?%k�>�r>�3�>�
|>D�>3Gx=�rS=��=O�/=�D�=�=R<�AD<��g;�/�;�Z�9�v9�,4�4��B                                        B���B�(#>+�>*%<�|�<�[%C��&C�-a?�?%�=�S=��2�[�2���)��y(�HdD*��D:�VA�A��?61?B�/X�N/X�N        1�ũ���{@��{@��                                                        {@��{@��B �<B, B �<B, {@��{@��@N�@T�88]�!�F���㡔��"�j��e�./�.�^.���0Z��6r�?,]�?+��-��-=Z@I�
@K<t.�K�/��m3_tg4�V&��񯯜C>��b?@[A��jA��A#d�A��?���?�e�=,n�=`��?���?�e�{@��{@��9}�>Re�9}�>Re�9}�$>Rg�A�e�A���A�e�A���{@��{@��C��4C���C��4C���{@��{@��C�	�C�	�?xJ_?|��?�  ?�  ?9�?S"~?Oz�?c��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�                                          ?U�?BcO?8�?S��?X9�?q��?l��?|�l?y�?}P�?z!�?}m�?zz�?}�=?{�?~�?{K�?~3t?{W�?~84                                        @,�E@,k?�u�?��N?��?���=��A=�~�C,��C6�6�T86�4�6�S�6�W@��;@��7>1&�>:dR?���?�{VB�ףB�٩=EV>���B�ףB�٩{@��{@��Bݗ-Bݗ-B+�B+�C[hC[hAHX�AHX�B��B��B�PGB�PGA���A���Bo�-Bo�-B�3B	y�@M�E@Z�iB�3B	y�B�3B	y�{@��{@��A���A��3@��@C��@��@C��{@��{@��>�I�>��j>�f[>�PA��+A��A�&A��A�^�A�o1@HU�@J��@8��@G�A/
cA;l�@�0@(�                                                                                                                                2�[�2���2��?3��32}�4 ��3��}3�^4�"3�R�3"�$2ɩ�2�<1n ./���/�:*w'�(+�z(m�"(��#&��'2�                                        0�!�0��,-��-WN[2Ӱ�2�7o2�ҡ3��36�Q4/�3�J�3�ۜ47?3���3#��2�H2]�1n�V0"`�0�\�*wN�(+�p(m��(��s&�5�'2�                                        N�c�N���B�p�B�zGE�.LE�p�8�D8�86�N6���5/��57D�2�W	2�=*��=0,�:�Z�@���A�!�A��?rLL?�>���?�>��?$�8>�a?8j�?[�?Q��?B�B?[�=?TH�?X��?Qq?R�z?R*�?R�?t+?zm>��$>���                                        @w��@���Fc�PFd�&Fc�%Fd��        7bz�7mW9;9J<8զI8塯8��*8�7�8��8�a7r��7`�6�ȑ6t��5?��5k~3�45�2�2O=R/_�?/���                                        A$�&A+��C�!;C��Z>�U�>�:0c��0^jq@m�Y@�[@I��@My�@و@�)?���?���>��>�d>�=н~<֎Q=�o;��S;���9��:ȿ75!'7�ۃ                                        ?�n�?��>���>��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��3��Ҵ=        C[�C+�6��B6��6���6�%I6�}�6��H@���@���>�W�>ߚ�6�V6��K        7]�7l{6�(6�X�A$��@��\6
&z6	��CTBC�+B�<�B��B]=A���A��@�7c@5��@�]<?�<	@A??7O=��=��i<��K<�;-:��:�d�                                        >y��>SJL�����3r�2��2oR�2l�01%}�0�1}/8t(.S�r,�D�-�,�&,�)+(�k+=Fu)Y�g)e��'sn'~��                                        3�*\3�8�@��@~m7@01�@
��?pEC?u�>V�>\�=���>��=$=�=b��<N�<c�;%��;.�x9˿�9�Z�7���7��3                                        <��<i&��з�غn1:6�1ɵ0��0�J/5+�.�x-B�9,��S+�+�*f*�?�(�n)jy'�'A*%#�%+�Q                                        1�a1Ȝg>��+>��>M��>!S�=���=&~�<U�x<`�;{(�;�)�:�k�;-E�:�k: �8�b8��7�i7���5�+�5�{%                                        5���5��uB��yB֡6S9�6b��D���Dw�.D<�KD0.C�HC��Cj�C9�
B��B�ynB'vIB�ȚA��lA���@�'�@���?V�?5<9�%<=��                                        ?�n�@@���%}�b�3gQM3�3Q�3"L72=��2:�\0��0p�s/+r/~z�.�4X/"w�-��%-�e�+ݿ+�)�}�)�?                                        3��4"�A���A���AbT+AU�~@��/@��_@Bg�@n�?���@Vk?U��?���>�sL>Ϳ�=��F=��3<A=_<Ha�9�u�9��:                                        =�8>
�=��ڲ� �1�֑1��1�c�1k �0e��0`q�.�	.�v,-,�-��@,[�[-��+t+� �)�/$)Ӎf'��P'��H                                        2q�2$9!?Ӭ{?�d�?���?��>��)>�<�>:֋>f4�=��	>c�=80e=��o<��<��D;���;�U�:
��:ʿ7?�7Cj�                                        6�,6�uCRC1�Y6[�6!�D��D���Dq|%DS|tD1D��C|OC��Cx�Cl��B�C �B2	-BXYA)s�A;��?J�?US<B�L<G̈                                        @:� @_;E�`�4� c�25�2(�383b2��2sێ2D�J1mO09Ԭ/*.��/ �l/�ac.*�%.Y��,���,���*8�*@5�                                        3zԓ3��A��LA���A�!�At�A!t�A"�Q@���@�PN@:��@�:�?�jD@U�e?k��?�!I>hH�>}mu<��}<�O99�(D9�I                                        >?�$>\����eH�~`0�r0ҕ�1�l!1=��0��W0�Ɵ/)��.�|�-b�|-W�s,�|z-�c�,��,5t�*�Ms*�0k(��(
��                                        1�$�1��?�"?£`?��?��?5�Y?6B>�>t>�p�>$��>��=�<�>1�==��=gV&<1c <D�M:Q��:]��7IN�7Nͱ                                        5wI�5��0C~��C�uI6�(6�X�7cTh7o^@f�@_=)y/==�H:4Ň:J��A �FAr�>$��>8�;/�/;E877?�X7H}�                                                        ����4=4�|�4�RW1zY;1zY;3��Ҵ=���i��w��.� 4"�M4�M54��2�N2�ڎ0�^803�I2�UE2ܤ07a��7o�|:�C:��<i�'<t	K<A�Y<<�;��;�v1;��;z�:�UM:���9�{�9�8�߷8��F75�7\�d5?X�5}4�2B��2���                                        ,���,���?�?xɎ?�o?e�k?�  ?O&?m]?.B?D�>�4P>��>V)�>`>OG�>Q��>N�v>M>L�>L��>L��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��CRxC��G�a\G�a\{@��{@��A۲A۲B�  B�                                                                                  B  B  Ap  Ap                                                          A�  A�  A�  A�                                  2�ev2�cS-�Ȩ-�Ȩ3S�2� o09��/�st7�T�7�x�4�W4���4�W4���2�}�2Σ2��2�~k2��2�~k5`�65+`1��1��N2�;2�32�G�2�V=��=��-k5�-k5�?04�?04�>��>��>_�>_�=We�=We�<M�l<M�l:��:��8�TB8�TB5B�5B�/���/���'`l�'`l�                                                ?"�?�7�I�7���?�Q�?��7�b�7�̩;�i�;�\�;�.;�.8N:j��6y=�6uM                7�;B7�%)        7��7���        ��s�6U�2            6�-�1\c�7`L�7��49�'6��W6�x�        6E4�    6E-�            6�T�7���6�T�7���6�T�7���{@��{@��6�36��K    B2�    B2�6���6�"5��5�k.74:�71� 7��;7��;@�YV@���3�`I3��:>�+{>я�1��2 B�VhB�#:B�VhB�#:{@��{@��7�7$39A�nkA���?��?�(�B�|B��o>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�>v�>:�@�@_/\L/P��+�(�+r��A�t�A�-�A��JA�mhA7j�A.L�@�a�@�;@"ڷ@9G?,��?��>2�>2E<��<Җ�:�o�;QO7�08:��                                                >�s?[�d4���4��I4���4��I3�<?3��94��4'�سNW.�M�43v�4?�~��N�D봌'ش�(1@2��@)��+s��+c��2�`�2�W�2�`�2�W�/����j�OE��G-�����2�/ݏ/���
 ����f3��\gA�.�A�ΚA��JA�$�AXk]ANȵ@�O@�FP@A&C@!M{?N5�?"�> �>7��<�;=�;��;-�N8,��8{                                        >��o?Y�e>�?��o>q��?�
?0�(?�M?�-�@)�?\�?�G
>�0j?1��>G�m>��=L�S=��e;�GE<��9�<9{��                                        <��5<��.Be�-|P�/��    90��:Оg9���;��;�rQ;�sE=O�Z<ЭR=�pc=��=O;��m:���;���0Y:.�.Q $.�&�,�ï-�                                        4}pk4��t3�u36�:2I�:NǠ8��9&�=/g=* X6pe6��4��+4�_06�zS6�zS=�H�=��^A�vA�L>y@�>��N=�V�=�:$6��	6��6ՠ�6�RA���A��`5���5�(�C]OCb8B�L�B�0&B���B���BzBU�(A�\�B
RIA:�pA��@�+�@�{`?��?�"`>#]�>)P                                                ?�S�@6@����N��0���;y2�21���1�.�2!0��0��l/Z�;/� <.�f/B��-�2�-��+�S,�*
��*z�                                        48;�47{A-ѿA5�JA��A$ @���@���@L��@��p?�&R@8m�?y�?�z7>��3>준=�!l=�-�<Y�"<as                                                =S��=��@�f�ܓ.[�q����0�1/+�d/l��/�8�.P�3.%
,�P- @,%b�,��++*�+IJ)��)�6�'��'��                                        1ă�1öx>�hC>���>�d>�,>F�\>��O=�G4> =i�z=ĹL=��=_qL<\�7<|k3;F�;U�@9�n%9�7�T7$                                        5�b5�C ��CR�L6'�6/25D^��D[,�DCY�DA�AD�lD�C�0�C�AC5F1C���B���C\C�Bt�PB�U�AmA�;�?k|�?wS�</<�                                        AV�^A������D��l׳��.�33��M3���3l�w3D�2 b17�^0zt.���0$-0��/m��/�_-�/�.
<�+E��+O�                                        4�"4�J�B�e�B��B�;�B�H�B0��BAZ=A�AQB ��Aq��A��A*w�A���@�,�@��?�`?�O�=���=��:H&�:MD                                        >���?育b'�Q�)�dT(��Sj1Jl1fL0�\;0Ѫ�/��5.���-���,	|-��.uA�,���-'��+�N�+�u�(Ҭ�(�WS                                        2\�2f� @J�@�g@
�e@	�K?���?�>?Mx?�T�? ��?[�>��?��>.�>X�B=(�(=<�;'|t;/��7��08 ��                                        6��6=�C��C��J6WK6�`DO��D'o�DGY;D!Z�D)�CD�NC��C��C�ȌC��CF��Ct&QB�*C��BPd�Bk�Am$�A��?�M/@$                                        A�nqA߫���������g�-�)jZ��Yy2H�4�LL4~�N4H3��2��.2��1�2uV�1d�1n��0#�0S!�.�P�/@3�                                        4��5�rB��B���B�z�B�BB��iBdIB<�zB3E9A�GB	�"A�ܛA�Q�A6j�AO�#@��@��?���?��>	�Y>{H9                                        ?�,?&C��вÏo��g˲X���hI/�e1س�1���1*��0��/�Զ/�)�//�I.@�K.��-Qį-��,�s,v                                        2#2+<�?Ԅ�?�t??�"?�9�?��?��C?q��?ewS?`�?0,�>�W�>�`>i~ >���=�e.=�<��o=#ŀ;0{;���                                        6/1�67�=C�+C��CeRCV�1Cd�CU�C_Q�CR��CSSdCJ�C?�/C>��C)�C.grChC*�C݄C
�B�m1B��B�u8B�[F                                        B:kVBH�µ$冀)O�?A��,���1�_ ܯ�c�2:LC3��3^��3bk�3"42�-�2�U1�o2�1�1��0���1 l�                                        4aG�4k�A�t�A��[A�v�A�A���A�w�A��A�U�A�;YA�zA��A���Ao�AxD%AQblA\��A5��AF�A�`A1�                                        ?n��?����SI�K���2<|�E����汎�Y���/nv$1�,0�y�0���0=�l/،�0	,6/�=/Am�.%E�.F:�."]�.$b                                        1�-�1�y�>��	>���>�	>��>��>ף�>�e�>���>�#>�+�>��>���>���>���>��>�nk>h�>~>E�>b�N                                        6��6�Dc(�D~a6��r7��6��r7��?�[�@=3p@8�@���@�f�@���A	fAd��Aq�JA�?�Ax��A���@�\B@�Q<                                                                @s�~@�d�@��	AN��Au�A���Bu�B7}B��B��SC
C�_C�s@C��
C�a�C�3�D��D�FD7�D7�                                        �a���h��`i��]���Al��Y���&�f�B�l��Hh��\ʿ+���Q�	��Q:�?�p�� ���~��9�����������p  �p  �p  �p  �p  �p  �p  �p  �p  �p  B&b�B��T:D�:'r;���;�uY;yP�;~0;"��; �:�p�:�k�:��9�l'9$�&9 Fr8�s8�Z6�fw6�^5
�J5<�`2^�k2�W                                        6��r7��        7��37��&D��
De��BOxB*�f@P�<@Wj        0e��0[E�?�?��{@��{@�μm���s�i        C�uC�u{@��{@��C��uC���C��;C��{@��{@��C�,�C�U%C�uC�u{@��{@��@z@]�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@���  �  ?�,?z4H�  �  �  �  EG�EE0pCp�C>NAC�AT�C~��C�uICx�C��d@�F�@��G>��>�ߔ>�PK>��]D�_�DȼB?>BQ�@~G�@���Dc(�D~aD��D�B���Bä�?�|�@��?�^h?�� D�KMD�3�D(1�D7�B?>BQ�@~G�@���C���C��C���C��{@��{@��C�1*C�)�C�1*C�)�{@��{@��C���C���?W�>?]��C���C���{@��{@��C�(pC�4�C�nC�.�C�
C�NC��C��C�ؚC���C���C�ШC���C��1C��fC��ZC���C���C��5C��TC��hC���C���C��ZC���C��6C���C��cC���C��|C�AC�${{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�nC�7E�hE�{@��{@��E�E�uE?\�?��?�?R~?�V?�6? ��?�g>��W>��d>�7�>��U>�C�>�uJ>�C>�)�>��>���>��C>���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��?|A�?~HA                                E�@ E�@         >�y�?)q�?b�?b�C��&C�c�6A[6-��6^f6-��                ;�r:��X?9?48?Q�?"��?9�?+�?�u?5k�?6�S?F��?Y��?\l4?t�??r�2?}W?~f?�c?�>?�J?�o{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��6o�5�\e�]u�n��A�  A�  >x��=��N>�5�>��?"?ī8{@��{@��{@��{@��        {@��{@��7l��7��|7B�l7Ly�7#
7 ��6��6�tA5��Q5���4���4|i3�y
3���2��2���1S��1`en/�۸/�[Z                                        {@��{@��