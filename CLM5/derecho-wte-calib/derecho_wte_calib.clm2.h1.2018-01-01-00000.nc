CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          (   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 04/16/24 10:35:17   source        Community Land Model CLM4.0    hostname      derecho    username      marielj    version       release-cesm2.1.5      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       derecho_wte_calib      Surface_dataset       Esurfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_c230220_v3.nc      Initial_conditions_dataset        .1ptBGC_spinup_grass.clm2.r.0201-01-01-00000.nc     #PFT_physiological_constants_dataset       #clm5_params_bayesopt_wte.c171117.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            (ltype_landice_multiple_elevation_classes            ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice               (ctype_landice_multiple_elevation_classes      4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      day_1      Time_constant_3Dvars_filename         /./derecho_wte_calib.clm2.h0.2011-02-01-00000.nc    Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE       "   levgrnd                	long_name         coordinate soil levels     units         m         d      �   levlak                 	long_name         coordinate lake levels     units         m         (       T   levdcmp                	long_name         coordinate soil levels     units         m                |   time               	long_name         time   units         days since 2011-01-01 00:00:00     calendar      noleap     bounds        time_bounds              �   mcdate                 	long_name         current date (YYYYMMDD)              �   mcsec                  	long_name         current seconds of current date    units         s                �   mdcur                  	long_name         current day (from base day)              �   mscur                  	long_name         current seconds of current day               �   nstep                  	long_name         	time step                �   time_bounds                   	long_name         history time interval endpoints              �   date_written                              �   time_written                              �   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��             �   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��             �   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��             �   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��             �   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����             �   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����             �   nbedrock                	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����             �   H2OSNO                     	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��             �   H2OSOI                        	long_name         0volumetric soil water (vegetated landunits only)   units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      P       �   QINFL                      	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !8   QOVER                      	long_name         surface runoff     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !<   QRUNOFF                    	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !@   QSNOMELT                   	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !D   QSOIL                      	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !H   QVEGT                      	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !L   RAIN                   	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !P   SNOW                   	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !T   SOILICE                       	long_name         #soil ice (vegetated landunits only)    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      P      !X   TSA                    	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            !�   TSOI                      	long_name         +soil temperature (vegetated landunits only)    units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      d      !�   ZWT                    	long_name         ,water table depth (vegetated landunits only)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            "   	ZWT_PERCH                      	long_name         4perched water table depth (vegetated landunits only)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            "<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  C�A�B>@�B��?�           E� 3�      	�     �@��     @��     04/16/24        10:35:17        B\�?P>���>ۢ�>�o�?��?� ?#��?5��? ��?@?f�?>�'�>�M>���>���>�`�>��$>�>�B���            4���        7�:A�Ana�A"�                                                                    Cr�pC���C��xC��3C���C�V�C��6C��8C�W�C�C�~FC���C���C�C�	~C��AC���C���C��~C�^C�5�C��C��C��`C��oC���A	��<#�
