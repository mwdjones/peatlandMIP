CDF      
      lndgrid       gridcell      landunit      column         pft    @   levgrnd       levurb        levlak     
   numrad        levsno        ltype      	   natpft        string_length         levdcmp       levtrc     
   hist_interval         time             title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 11/27/23 15:59:35   source        Community Land Model CLM4.0    hostname      ubuntu     username      mwjones    version       5823c39    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
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
_FillValue        {@��   missing_value         {@��           �;�r<���=�=�o�>YI:>�l�?�~?��?�#'@7U�@��,@��rAN�A���B��=L��?��@ff@�33A��AI��A���A���B	L�B3�;�r<���=�=�o�>YI:>�l�?�~?��?�#'@7U�@��,@��rAN�A���B��º�Bº�BB>@�B>@�EA4�EA4�        ?�  ?�              FG�  ~�      1�     ��@��     @��    11/27/23        15:59:35        4���4���1��1B�/�F1.Wp�6��*6�U`5���5��B��B��B��B��B��B��7���7���7���7���1��1�:<                *�ۥ*�ۥ6�6�6�B�        ?Dv�?B8�        3�Y4�qf3e�3_��1���2�lT0�R�/�{�/��k1�4�1���0���    .��W        =���<���>G�=<�>6vj=6�b>8=So�=���=�Qj=u�=|��<oT%<�	;L�;��9���:2�7�(7���                                        8EK�;���8_z	<BAe8�!<���9�gI=�2<O�=��=/�6=��<<��=p;�SD;���9��,:>17��7�TN                                        >`ͭ>3R�=���=ϼ�=��3=�$2=�@�<��<�:޽�;[No4�<_:���7�_{9GP	4�~�0O>X+�R+�+�                                        @��@BV�@��@�M@��4?�#Q@2�?A0?�l�>W�=�}�:��:���4t_�2���,ocA+/X+!�+�+�                                        D7˗D8J�D ��C�hd        6��5��H5�v�5�5�4�Z#E��E�}/E�5�Eu�E�GE	Dn2�D��&C��C�i�CZ-C�eA�A�A�]@5�?@.��=��G=�w                                                ?�[�?��p1Z
�1H��0���0}�0Av��AU��A��A)�@��@��?�u�@��?cA
?�J�>���>�6�=�n,=~X4;���;���9�*9��                                                >0�n>.r0a$0Į/3�/%.U@!#�@O?��E?�[+?8��?;�\>�I>���>�i>)6=P1"=X�T<(�S<$��:w��:n��7��o7�c�3�',3��B                                        B�ءB��Q>9�>�<�1�<�gC�C��?W?	�=�;r=��W2��0��k*�~�)�0�D��D�@�7@�"L?��?�/X�N/X�N        1��O1�x{@��{@��                                                        {@��{@��B��B+{B��B+{{@��{@��?�"�?�8,7��(7�,%?9ȧ A%����"��"�+�2$�.���.����Cȥh�?1�?\�-�%i-`u@78�@(}�.��^/Ec3Pw�4��䪔�D�D��>�n?�-AĆ0A�ԢA@EAR�F?�*A?ߤ}=�_�=�BC?�*A?ߤ}{@��{@��9��>>d��9��>>d��9��>>d��A��A�{�A��A�{�{@��{@��C�jC�a�C�jC�a�{@��{@��C�	�C�	�?n��?k�7?�  ?�  ?$S�?#�?8� ?2(?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�                                          ?
v?Y?#�o?%r�?L��?>?bK?Qܗ?g��?b3?m��?k�%?s?q�?s�	?qʤ?s��?rU�?t�<?s~6                                        @ �@��?���?��?���?��=���=_��Cw`C6�Z�6�̓6�C6���@`��@J�m>h>K�?� ?�+RB��0B��C=(�>�R�B��0B��C{@��{@��Bݗ-Bݗ-B+�B+�C[hC[hAHX�AHX�B��B��B�PGB�PGA���A���Bo�-Bo�-B	��B[z@u��@�ZyB	��B[zB	��B[z{@��{@��A���A� V@�`@FS�@�`@FS�{@��{@��>�j�>�mi>���>�A�/�A��A[�A�(A��A��$@HG@J�@:2�@K�3A0�A@��@!��@-.z                                                                                                                                2��0��k0�j�0ؗ�1Ў1Y9n3�%1��3��2� 2a�1��1q:�0��/0�".>�)��`)N��%r4$��!{{�!\V�                                        /��k.{d^,̤�+s:&2&��0��R1>F1ʚ1�g�1tqZ3 �1���3���2�2c��1��1��D0�T/]�.��)��)N��%��$���!{��!\�%                                        N�*�N�,B�B���E�(�E�N�7��Z7��|6^(�6@�5
B�4�� 2��u2h�=`Z<�ӄ:���@��vA�y�A���?k�?�(@>���?#oj>��?.�8>���??��?�0?S��?Dt*?\Z??UH?Xp�?Q�?RD�?R+i?R�l?tJ?yF>��d>��$                                        @y�t@�EFc�tFd��Fcz�FdyW        76��7 ��8ړ�8��8��8��8f�R8E�7���7�پ77�7�6&%V6�4�Ϩ4��2��)2�X^/[k/H��*E***O�                                        A!��Ac�C�xC��>ʬ*>ƫA0RkN0NCO@]��@x��@A��@L7�@�?�7?��?�>�&>���=��=�1�<9��<0�:c:��6���6�&02 �1�u                                        ?�Ï?��0>�}+>s�.{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@�ε��ε�I�        B��B�}�6`\Y6FC>6a�a6H�m6\d 6C�'@b�@N$4>��>���6s��6`XW        7:�;7&�]6�Pv6��1A.tA�5�X5��C?ۉC��B�q�B�<�B��B��A$AIO�@_� @�,m?�w�?��C>��>��z=�e=�Ic<x�.<q;�9V�9I+�                                        >g`<>V@L���5�t
�2�Ӆ2�c\2J��2��1�00�6+/Hm�/d?-H=K. O,�^!,�r�+	+�z)2u�))�V'j��'c�#                                        3�v�3vx@��@Z/6@X:@�C?hY�?\&	>uU>���=�ߎ>�9='��=7�<)xp<'�;{�;\9�ʣ9�L6�H6�i                                        <}ϭ<h�ɱ�4I��nS1*�0�;�0|��0F�y/�v/
��-M�?-&��+Xs�,Y�*�U?*��c(��(��m&�&��z%o%Xt                                        1�[�1��>�+C>�ƿ>5T�>((=r�f=jE<k�<�c;��{;���:�F�;|�99�$�8���8�@|7{��7tB=4av(4S                                        5���5��B�Q�B�F6-Y�6.�D�~DW�aD/k�D�"C���C�NC�uC�B}�B�uBGB7��A8_�AL�a@(j@(�>�.>�:3�:-��                                        ?ˆ�?ԞĴ��\�-a�3O�2�`�3:(72�K�2(�2 �40�D�0��]/n�/L��.�¦.��C-G�!-k��+��+���) ��(��n                                        3��3�z�A�W'Ak��AH�oA00Q@��S@��@-�@I<r?�~?�B??A�?���>��O>�V�=}�=~a�;KD8;E�!7�Pm7�BT                                        =���=��O�ͿV����1��1{h1��>18$F0PR0G��.�v;.���-6�B-u�,x/,��+l�+,?�)�/))���&��7&���                                        1���1��A?��?���?{|?[5{>��->��>(��>C!=��=ɰj=��=Ig<A� <W	e;-}
;-�19	9le5:\754��                                        6S�x6=�B�E�CU�5�\>5��wD�CFDxׂDgz�DCKC�+C�_�C]��C�\|B�yC1�OB�7�B��LA�@fA��@M��@J0�=��a=�59�S9�G�                                        @J�@,�0�T}���s1���2�3<�	2�q@2c��22��0訔0X,.�yH/�/�/�ة.��.)G�,	�+��h(���(� P                                        3CG&3.�A�?hA�(A���AW A��Aݠ@�<�@��]@%��@s��?��@n!?�?.k�=�Ύ=��;/&:��|7�d7�                                        >#�)>+�V����oE�0��I0��1��13��0�7L0u(/Q).��=-_��-�-!w-f��+�[+�)�Y�)�a�&b)R&Z�R                                        1^+l1H(�?��G?��?�)s?�O?%B'?!�>��C>�� >0�>O/Y=�!5=��<��h<�
�;T�I;Q�8�~�8���4�z�4�x                                        5C��50[AC\>�Cf�j6�Pv6��17;#7'v�?�^�?ɃJ=Q%= ��:�:	��@�b2@�р>	4�=�^;Z?;��7&�I7l�                                                        5���5�I�4��p4n�1zY;1zY;���ε�IԳ�Յ���	5�S!5��4A�#4��2\z 2=��0��)0J�2��2�}�7P��7:9v8�8v9�D9�e9���9��96�89#��8�N8��t7���7�8�6��6]4��4��
1��1͢.�-��^)�<.)�B�                                        ,���,���?�?v�d?��?b��?�  ?E��?:J?%5?<X�>��'>��&>`�>p >Py�>WϺ>O��>M<�>L�e>L��>L��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��B��G�a\G�a\{@��{@��A۲A۲B�  B�                                                                                  B  B  Ap  Ap                                                          A�  A�  A�  A�                                  2��|2��-�Ȩ-�Ȩ3;�2>*
0$/<�7p�97P�a4�/C4�r4�/C4�r2�}�2��\2���2��,2���2��,5c/5@�1��y1�#72^��1e�2(:�0��s<i�<i�,�=�,�=�>��$>��$>Ja>Ja=�x�=�x�<��H<��H;���;���:/�:/�8�8�4��a4��a/H��/H��&�`�&�`�                                                ?* ??8N7s��7V�?į�?�YV7�
u7w��;���;�X�;�.;�.8�A:s�E6�0;6�
                7�1�7��D        7�/�7{�        �ف6h�            6���0�c7n!@7��P9څm6��6|2R        6E)�    6E7�            6��C7�S�6��C7�S�6��C7�S�{@��{@��6�<`7�S�N/p_�N/p_6���6�ޛ5���5�H�7$�_7p<7��;7��;@M(�@9�93�_3p��>��q>��1��U1��)B��B�8FB��B�8F{@��{@��7��6�thA��A�^I?�O�?�FB���B�F�>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��=�Ɣ>��>���-���-ތ�*
�(*x2@^[�@b�@)`�@%��?Υ?��;?U1?>��>�6b>u%�=o~"=''8;���;��99 �9�5u�5L�0�ܑ0��                                                >�=�=�4WV42A�4WV42A�3���3��4�3�Wȳ0(a�V�4�54 ���|���^��h��Ks?7_Y?/�+l �+]�@2c��2=J�2c��2=J�0\0�N��s}�ɵ����ϰ��R/�{�/�@,�ߡ���9S�����52@�!@�D@��(@�s@h�6@a"j?���?�v?9G�?B>�&=��<��a<zh :6^�:,�*6Ɗ�6���2s�2�i                                        >�=��=,W�=�y�=xGh=g?[>��w=��	?,_7>p��>��>��=�2=stv<�(R<�fi;9�);��8*)8;�3���3`[\                                        <<(�:^��-�0�+�Pv/��    7�TJ7�M�8<-7�;<-|8�-�=nH�;���<�
�;=G�<|'9I��9EZ�9Q��/~��/5��+��*���']�'A�f                                        4t��4��3�P38�*:/q�:H&�8��9�=/��=&��6 Q�6OU4��4�f6�zS6�zS=��R=�*gA}��A���>o��>��=�XH=ï	6���6���6��G7�aAu\�A���5���5��OB��#B�;uB��sB���Bj��B���B #B�wA�x�A��A ґA[��@e
�@|��?Wu3?V�=}�=@                                                ?��M?�/���|�8B:1�pY���2E��1��z1���1���0�/�0��/%�/1�N.�
�/w�-a�]-�W�+磇+��(��W(�$                                        4H�4_�AMmA|�@��@�f@�o�@��2@*�z@K��?�K�?���?Vn?�Y�>���>�`�=��w=�Hu;A��;<Uw                                                =.{�=E��*�@�Ċ�/U�,����/ҼW/��/K��/E �.&�.�,�1,�թ,.��,�n�*�+�)w�)v�&���&kQ                                        1���1���>��A>��>{uC>r��>&ݤ>7�=�=O=��=G�t=���<乱=�<"��<3�/;6�;��8�T�8�E�4���4��,                                        5I4�{B�#TC	_6
�w5�?�D9E�D7�D#D�D�vC۞C�Ck�C�pC�|CV�BŒ&C�BDB�!@P�u@I_�=t|<�HX                                                A �8A6�*��1&���T���S���i3��w3Hu63J-03	A1�/#1#%�/��/��0=�;0���/:�/S�,�,�/(���(�`r                                        4���4��Bw�BTJRBY�>B;��BieB�jA���A���AD�A�P�A��A=@:Ѱ@Nl,>�N>�?�;0�P;&��                                                >�y�>�v�6�Ӳf��[�&����1.w�0��\0ק�0�9�/rTj.��-|_�-�;:-�1Q.9/�,�c�,�+~*|rC*s��&px�&cB5                                        26��2%n�@��?�qh?�4�?�S%?�,'?��?'$�?G/�>��?ޗ>�~�>ɝ�=�F=�/<M2<3A8�KF8���4;Q�41��                                        5�P5�	DGZD'2�5�==5��E:�gE06[E44E+�E>E�D׵�D�iD�&4D��pC�D�%B�s)B���@kRk@Nf�=��<�5b                                                B}>�B��?�jҪ�:�w�I	Ү� ������4�u-4��3��=3_v�3y83Eؚ2W*�2u�0���0�4�-���-���)���)|�                                        4�Ĩ4��JC�N�C��HC�)�C��CxcQCt�:C,��C=] B�ֹB�MB@�5BP�A
A;>�A�>��;oo�;J�                                                ?��?�6��Ib�n�p��<}�0iձ%���m1��x1&�0�_0�z02��0}>6/���/�S6-�P-�v�*ȉ*���&�Q�&��                                        2�1��@�@�q@��D@�5�@���@���@\�L@rb�@K�@ ��?vy�?��.>1��>&��;��;�[�8�+"8��z4b^�4<�c                                        6�45�O�C�K�C��gD:$[D6��D8ٕD5MD3�~D0l�D&#kD#��D�FDM�C��C�\C\M
CU!/B�IfBw�@�W6@��h>�av>��O                                        B�B��%=M�&+x���Դ� ���<W��6���<E3A��2��3v[�3k>3J13VFB2��l2��1���1q�l/SuK/:�                                        4Q`3�_�B���B��B��DB�
rB���B�#�B��"B�1BB_�qB^�B"7KB!kA�=�A���@���@�E�?I�^?1j�<�5�<�%L                                        ?E�}?AZ��S���T���9�%�,bW��{����Z��R���0x:�/Ћ0��40�m�0�Z]0�"�/�/���.��..��",�T�,n6                                        1Fΐ1Q�?��?���?�IL?���?�	?���?� ,?��"?�.]?���?O��?N1->ᖒ>�>�>	�w=���<�2�<cZ: 5�:�@                                        5�%�5���D�s�D�]6��6���6��6���?��@[��@=+�@�{n@���A	��AǪAn/6A��lA���A��A��Q@��A�    ;E8C                                                        @xF@� @쒯AS��Aw[DA��B�)B6�XB���B��C�C+�C�> C�8gC�b�C�
�D��D�?D7D7                                        �V���tbD�d\1�`Y��Y���e�:�*�C����-���`��/�R�p���m��凸��\���黭8������*��$�p  �p  �p  �p  �p  �p  �p  �p  �p  �p  B)�#B���9��9�/�;�%;��;o�;|��;5y;�5:�a
:�Ep:�_9��8���8�6=7di�7YR5.�65'��1���1���-]�-�                                        6��6���        7��7�hZD���D�4�A�ZDA�b@ <�@�}        0�k0K��?�0?��{@��{@�μkcԼk�T        C�uC�u{@��{@��C�2C�%C��YC���{@��{@��C�4�C�jQC�uC�u{@��{@��?�i�?��/{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@���  �  >��>�6��  �  �  �  Eg�EEdegC��CfN@�=�@���C\>�Cf�jCY�CdK1@��"@�6K>���>�{�>���>�7aD�3�D�șBQ	B�{@G�@2)D�s�D�]D��D��'B◘B��@��@!hD@�@?D�DD�;KD�`D�*BQ	B�{@G�@2)C���C��RC���C��R{@��{@��C�1vC�)�C�1vC�)�{@��{@��C��YC���?Bv�?7�]C��YC���{@��{@��C�0tC�IC�#�C�A�C�4C�/�C��mC�.C���C���C��MC��C���C���C���C��IC���C��gC��TC��FC���C���C��tC��
C��bC���C��9C��C��C���C��C�6�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�zC�*\E��E�'{@��{@��E�͌E�?3d?� ?��?�?U\?�<? �?�'>�G9? 3�>��
>�)�>�A�>�FA>���>�w�>�E�>��Q>��>�C�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��?(?�_�                                E�@ E�@         >��?5�?b�?b�C��!C��}5�["5ޤ�6,6
�4                ;#��:�0?A?X�?9�?$7b?!�?+�?m�?4a,?5��?B�'?W��?W��?tO{?p��?}�?}۟?�?��?�g?�Q{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��6���6�8�;v��.V*A�  A�  >m7N=��>>�w >��?�k�?�k2{@��{@��{@��{@��        {@��{@��7E�e7?��7$=�7�6��]6�j86�H'68c05�5j8�4��4_�3o�3(e2kX�2/t�1 1�0��~/u��/A�                                        {@��{@��