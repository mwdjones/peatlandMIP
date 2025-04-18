CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          '   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 04/16/24 10:38:20   source        #Community Terrestrial Systems Model    hostname      derecho    username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        hillslope-wte-calib-pfts   case_id       hillslope-wte-calib-pfts   Surface_dataset       hsurfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_HAND_3_col_hillslope_wte_pft_calib_fmax0.1006.nc   Initial_conditions_dataset        finidat_interp_dest.nc     #PFT_physiological_constants_dataset       clm50_params.c240105.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         -./hillslope-wte-calib-pfts.clm2.h0.2011-01.nc      Time_constant_3Dvars      AZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE:PCT_SAND:PCT_CLAY         W   levgrnd                	long_name         coordinate ground levels   units         m         d      AL   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m         P      A�   levlak        	         	long_name         coordinate lake levels     units         m         (      B    levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m               B(   hillslope_distance                 	long_name         hillslope column distance      units         m               B,   hillslope_width                	long_name         hillslope column width     units         m               BD   hillslope_area                 	long_name         hillslope column area      units         m               B\   hillslope_elev                 	long_name         hillslope column elevation     units         m               Bt   hillslope_slope                	long_name         hillslope column slope     units         m               B�   hillslope_aspect               	long_name         hillslope column aspect    units         m               B�   hillslope_index                	long_name         hillslope index             B�   hillslope_cold                 	long_name         hillslope downhill column index             B�   hillslope_colu                 	long_name         hillslope uphill column index               B�   time               	long_name         time   units         days since 2011-01-01 00:00:00     calendar      noleap     bounds        time_bounds             D�   mcdate                 	long_name         current date (YYYYMMDD)             D�   mcsec                  	long_name         current seconds of current date    units         s               D�   mdcur                  	long_name         current day (from base day)             D�   mscur                  	long_name         current seconds of current day              D�   nstep                  	long_name         	time step               D�   time_bounds                   	long_name         history time interval endpoints             D�   date_written                             E   time_written                             E   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��            B�   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��            B�   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��            B�   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��            B�   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����            B�   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����            B�   nbedrock                	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����            B�   
grid1d_lon                 	long_name         gridcell longitude     units         degrees_east   
_FillValue        Gh��y �            B�   
grid1d_lat                 	long_name         gridcell latitude      units         degrees_north      
_FillValue        Gh��y �            C   
grid1d_ixy                 	long_name         ,2d longitude index of corresponding gridcell   
_FillValue        ����            C   
grid1d_jxy                 	long_name         +2d latitude index of corresponding gridcell    
_FillValue        ����            C   
land1d_lon                 	long_name         landunit longitude     units         degrees_east   
_FillValue        Gh��y �            C   
land1d_lat                 	long_name         landunit latitude      units         degrees_north      
_FillValue        Gh��y �            C   
land1d_ixy                 	long_name         ,2d longitude index of corresponding landunit   
_FillValue        ����            C$   
land1d_jxy                 	long_name         +2d latitude index of corresponding landunit    
_FillValue        ����            C(   	land1d_gi                  	long_name         '1d grid index of corresponding landunit    
_FillValue        ����            C,   land1d_wtgcell                 	long_name         2landunit weight relative to corresponding gridcell     
_FillValue        Gh��y �            C0   land1d_ityplunit               	long_name         Clandunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)    
_FillValue        ����            C8   land1d_active                  	long_name         (true => do computations on this landunit   
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          C<   
cols1d_lon                 	long_name         column longitude   units         degrees_east   
_FillValue        Gh��y �            C@   
cols1d_lat                 	long_name         column latitude    units         degrees_north      
_FillValue        Gh��y �            CX   
cols1d_ixy                 	long_name         *2d longitude index of corresponding column     
_FillValue        ����            Cp   
cols1d_jxy                 	long_name         )2d latitude index of corresponding column      
_FillValue        ����            C|   	cols1d_gi                  	long_name         %1d grid index of corresponding column      
_FillValue        ����            C�   	cols1d_li                  	long_name         )1d landunit index of corresponding column      
_FillValue        ����            C�   cols1d_wtgcell                 	long_name         0column weight relative to corresponding gridcell   
_FillValue        Gh��y �            C�   cols1d_wtlunit                 	long_name         0column weight relative to corresponding landunit   
_FillValue        Gh��y �            C�   cols1d_itype_col               	long_name         #column type (see global attributes)    
_FillValue        ����            C�   cols1d_itype_lunit                 	long_name         Jcolumn landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)     
_FillValue        ����            C�   cols1d_active                  	long_name         &true => do computations on this column     
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          C�   cols1d_nbedrock                	long_name         column bedrock depth index     
_FillValue        ����            C�   
pfts1d_lon                 	long_name         pft longitude      units         degrees_east   
_FillValue        Gh��y �            D    
pfts1d_lat                 	long_name         pft latitude   units         degrees_north      
_FillValue        Gh��y �            D   
pfts1d_ixy                 	long_name         '2d longitude index of corresponding pft    
_FillValue        ����            D0   
pfts1d_jxy                 	long_name         &2d latitude index of corresponding pft     
_FillValue        ����            D<   	pfts1d_gi                  	long_name         "1d grid index of corresponding pft     
_FillValue        ����            DH   	pfts1d_li                  	long_name         &1d landunit index of corresponding pft     
_FillValue        ����            DT   	pfts1d_ci                  	long_name         $1d column index of corresponding pft   
_FillValue        ����            D`   pfts1d_wtgcell                 	long_name         -pft weight relative to corresponding gridcell      
_FillValue        Gh��y �            Dl   pfts1d_wtlunit                 	long_name         -pft weight relative to corresponding landunit      
_FillValue        Gh��y �            D�   pfts1d_wtcol               	long_name         +pft weight relative to corresponding column    
_FillValue        Gh��y �            D�   pfts1d_itype_veg               	long_name         pft vegetation type    
_FillValue        ����            D�   pfts1d_itype_col               	long_name         'pft column type (see global attributes)    
_FillValue        ����            D�   pfts1d_itype_lunit                 	long_name         Gpft landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)    
_FillValue        ����            D�   pfts1d_active                  	long_name         #true => do computations on this pft    
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          D�   FSAT                  	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E,   FSNO                  	long_name         "fraction of ground covered by snow     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E8   H2OSNO                    	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            ED   QICE                  	long_name         ice growth/melt    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            EP   QINFL                     	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E\   QOVER                     	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            Eh   QRUNOFF                   	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            Et   QSNOMELT                  	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   QSOIL                     	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   QVEGT                     	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   RAIN                  	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   SNOW                  	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   TREFMNAV                  	long_name         (daily minimum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   TREFMXAV                  	long_name         (daily maximum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   TSA                   	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   ZWT                   	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   	ZWT_PERCH                     	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            E�   TSOI                     	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��     ,      E�   H2OSOI                       	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      G$   SOILICE                      	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      H<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  @      @]�;G�a�@j�7��@d      @�V��%�
@��	� @��     @�     @�T             ?�o��b��@�s���>�MP��>��I�+Ј>ƍH�s�@2ƂH@.���f@�l�:~         ����         ��������C�A�B>@�=xČ?�           @p�1&�x�@G�bM��      @p�1&�x�@G�bM��         ?�            @p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��                                    ?�3��ԣz?ڡK�6�?�h�i�2?�3��ԣz?ڡK�6�?�h�i�2                                    @p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��                                             ?�3��ԣz?ڡK�6�?�h�i�2?�3��ԣz?ڡK�6�?�h�i�2?�      ?�      ?�         
                                 D�� 3zu      �     ?P@�      @��     04/16/24        10:38:20        =�	�=� o=�[�                        {@��{@��{@��>IB!<��=8���8cG�6)�6*9�7�]d66�H6�Ny            7�
26LĂ6�,r6���7���7�.7�F7�;�7��6            C�3�C�bC��C�9�C�	>C���C�yC�D`C�@_<�u>�0�>��]@�R@�R@�RC��C�'+C�L�C�V�C�kC�'VC�
C��1C���C���C��}C��$C�'�C�*�C�N3C�f�C��"C���C�zC���C��C��SC�ҡC���C��C�ԜC��C��C��C�»C��C�t�C��C��C�;�C�ˤC��cC�&fC��0C���C�#�C��xC��C�0aC��EC�>�C�K5C���C�g#C�o�C��C���C���C���C��GC��aC���C��C���C�UC��\C���C�/+C��C� hC�E�C�8BC�B�C�XC�]TC�aZC�i=C�pYC�q7C�s1?m�h?m�h?�j?l1?l1?<Sp?ix�?ix�?iBI?e�T?e�T?e�T?`Ĝ?`Ĝ?`Ĝ?Y��?Y��?Y��?Tz�?Tz�?Tz�?Tz�?Tz�?Tz�?Tz�?Tz�?Tz�?Tz�?Tz�?Tz�?(�k?(�k?(�k?(�k?(�k?(�k>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��>��                                                                                                                                                                                                                                                