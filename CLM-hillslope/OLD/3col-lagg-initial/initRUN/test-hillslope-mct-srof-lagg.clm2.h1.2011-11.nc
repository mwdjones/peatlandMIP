CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          '   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 08/22/23 11:52:42   source        #Community Terrestrial Systems Model    hostname      cheyenne   username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       test-hillslope-mct-srof-lagg   Surface_dataset       Tsurfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_HAND_3_col_hillslope_lagg.nc   Initial_conditions_dataset        finidat_interp_dest.nc     #PFT_physiological_constants_dataset       clm50_params.c211112.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         1./test-hillslope-mct-srof-lagg.clm2.h0.2011-01.nc      Time_constant_3Dvars      AZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE:PCT_SAND:PCT_CLAY         U   levgrnd                	long_name         coordinate ground levels   units         m         d      ?<   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m         P      ?�   levlak        	         	long_name         coordinate lake levels     units         m         (      ?�   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m               @   hslp_distance                  	long_name         hillslope column distance      units         m               @   
hslp_width                 	long_name         hillslope column width     units         m               @4   	hslp_area                  	long_name         hillslope column area      units         m               @L   	hslp_elev                  	long_name         hillslope column elevation     units         m               @d   
hslp_slope                 	long_name         hillslope column slope     units         m               @|   hslp_aspect                	long_name         hillslope column aspect    units         m               @�   
hslp_index                 	long_name         hillslope index             @�   	hslp_cold                  	long_name         hillslope downhill column index             @�   	hslp_colu                  	long_name         hillslope uphill column index               @�   time               	long_name         time   units         days since 2011-01-01 00:00:00     calendar      noleap     bounds        time_bounds             B�   mcdate                 	long_name         current date (YYYYMMDD)             B�   mcsec                  	long_name         current seconds of current date    units         s               B�   mdcur                  	long_name         current day (from base day)             B�   mscur                  	long_name         current seconds of current day              B�   nstep                  	long_name         	time step               B�   time_bounds                   	long_name         history time interval endpoints             B�   date_written                             B�   time_written                             C   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��            @�   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��            @�   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��            @�   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��            @�   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����            @�   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����            @�   nbedrock                	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����            @�   
grid1d_lon                 	long_name         gridcell longitude     units         degrees_east   
_FillValue        Gh��y �            @�   
grid1d_lat                 	long_name         gridcell latitude      units         degrees_north      
_FillValue        Gh��y �            @�   
grid1d_ixy                 	long_name         ,2d longitude index of corresponding gridcell   
_FillValue        ����            @�   
grid1d_jxy                 	long_name         +2d latitude index of corresponding gridcell    
_FillValue        ����            A    
land1d_lon                 	long_name         landunit longitude     units         degrees_east   
_FillValue        Gh��y �            A   
land1d_lat                 	long_name         landunit latitude      units         degrees_north      
_FillValue        Gh��y �            A   
land1d_ixy                 	long_name         ,2d longitude index of corresponding landunit   
_FillValue        ����            A   
land1d_jxy                 	long_name         +2d latitude index of corresponding landunit    
_FillValue        ����            A   	land1d_gi                  	long_name         '1d grid index of corresponding landunit    
_FillValue        ����            A   land1d_wtgcell                 	long_name         2landunit weight relative to corresponding gridcell     
_FillValue        Gh��y �            A    land1d_ityplunit               	long_name         Clandunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)    
_FillValue        ����            A(   land1d_active                  	long_name         (true => do computations on this landunit   
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          A,   
cols1d_lon                 	long_name         column longitude   units         degrees_east   
_FillValue        Gh��y �            A0   
cols1d_lat                 	long_name         column latitude    units         degrees_north      
_FillValue        Gh��y �            AH   
cols1d_ixy                 	long_name         *2d longitude index of corresponding column     
_FillValue        ����            A`   
cols1d_jxy                 	long_name         )2d latitude index of corresponding column      
_FillValue        ����            Al   	cols1d_gi                  	long_name         %1d grid index of corresponding column      
_FillValue        ����            Ax   	cols1d_li                  	long_name         )1d landunit index of corresponding column      
_FillValue        ����            A�   cols1d_wtgcell                 	long_name         0column weight relative to corresponding gridcell   
_FillValue        Gh��y �            A�   cols1d_wtlunit                 	long_name         0column weight relative to corresponding landunit   
_FillValue        Gh��y �            A�   cols1d_itype_col               	long_name         #column type (see global attributes)    
_FillValue        ����            A�   cols1d_itype_lunit                 	long_name         Jcolumn landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)     
_FillValue        ����            A�   cols1d_active                  	long_name         &true => do computations on this column     
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          A�   cols1d_nbedrock                	long_name         column bedrock depth index     
_FillValue        ����            A�   
pfts1d_lon                 	long_name         pft longitude      units         degrees_east   
_FillValue        Gh��y �            A�   
pfts1d_lat                 	long_name         pft latitude   units         degrees_north      
_FillValue        Gh��y �            B   
pfts1d_ixy                 	long_name         '2d longitude index of corresponding pft    
_FillValue        ����            B    
pfts1d_jxy                 	long_name         &2d latitude index of corresponding pft     
_FillValue        ����            B,   	pfts1d_gi                  	long_name         "1d grid index of corresponding pft     
_FillValue        ����            B8   	pfts1d_li                  	long_name         &1d landunit index of corresponding pft     
_FillValue        ����            BD   	pfts1d_ci                  	long_name         $1d column index of corresponding pft   
_FillValue        ����            BP   pfts1d_wtgcell                 	long_name         -pft weight relative to corresponding gridcell      
_FillValue        Gh��y �            B\   pfts1d_wtlunit                 	long_name         -pft weight relative to corresponding landunit      
_FillValue        Gh��y �            Bt   pfts1d_wtcol               	long_name         +pft weight relative to corresponding column    
_FillValue        Gh��y �            B�   pfts1d_itype_veg               	long_name         pft vegetation type    
_FillValue        ����            B�   pfts1d_itype_col               	long_name         'pft column type (see global attributes)    
_FillValue        ����            B�   pfts1d_itype_lunit                 	long_name         Gpft landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)    
_FillValue        ����            B�   pfts1d_active                  	long_name         #true => do computations on this pft    
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          B�   FSAT                  	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            C   FSNO                  	long_name         "fraction of ground covered by snow     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            C(   H2OSNO                    	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            C4   H2OSOI                       	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      C@   QICE                  	long_name         ice growth/melt    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            D0   QINFL                     	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            D<   QOVER                     	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            DH   QRUNOFF                   	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            DT   QSNOMELT                  	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            D`   QSOIL                     	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            Dl   QVEGT                     	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            Dx   RAIN                  	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            D�   SNOW                  	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            D�   SOILICE                      	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      D�   TSA                   	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   TSOI                     	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��     ,      E�   ZWT                   	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            F�   	ZWT_PERCH                     	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            F�<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  @      @]�;G�a�@j�7��@d      @�V��%�
@��	� @��     @�     @�T             ?�o��b��@�s���>�MP��>��I�+Ј>ƍH�s�@2ƂH@.���f@�l�:~         ����         ��������C�A�B>@�B��?�           @p�1&�x�@G�bM��      @p�1&�x�@G�bM��         ?�            @p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��                                    ?�3��ԣz?ڡK�6�?�h�i�2?�3��ԣz?ڡK�6�?�h�i�2                              
      @p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��                                             ?�3��ԣz?ڡK�6�?�h�i�2?�3��ԣz?ڡK�6�?�h�i�2?�      ?�      ?�                                          C�  2�a      N      >�@s      @t�     08/22/23        11:52:42        >��=�th=��>�e�>;?�>A�4@d/@��@
��>���>�d�>���>��E>�P�>�T;>�6�>��)>�F�>��_>���>�j�?��?��?�o? ��?E�?n�?@��?$2�?$xw?Tz�?#�n?#�?Tz�?,c�?,�t?Tz�?�m?��?�?�`?�?Z??��?�|>�Cp>��>�ќ>�k>�2�>��>��>���>��h>��>�
>΃>�T�>��>��>��>��>��>��>��>��>��>��>��{@��{@��{@��5���6� 6�?15�=58	55�X5��:5.�156,�6
��6s˴6tN^5�5�J�5�d�N��I�/6[3�6[i6YO6��=6���6�W�@]YZ@O��@P��@��mA(�A(d�@~�>A�eA	-M    @�B@���                                                                                                                                                                                                C���C���C���C��tC��C�$C�
CC�a�C�`�C�s,C��7C��FC��	C�P~C�O�C���C��C���C�J,C�o@C�nC���C�}C� C��\C���C��C��nC��fC���C��C��/C���C���C���C��C���C�дC�ϫC��dC���C��C�g�C��C��EC�HYC���C��KC�$C�j�C�jVC���C�]�C�];C��.C�[�C�[�C���C�d�C�dvC���C�t�C�t�C���C���C��C��%C���C��C���C��GC��FC��fC��9C��9C���C���C���?[H@���@��?��B?g$�?e�	