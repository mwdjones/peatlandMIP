CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          '   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 04/18/24 16:49:43   source        #Community Terrestrial Systems Model    hostname      derecho    username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        hillslope-soiltest-pfttest     case_id       hillslope-soiltest-pfttest     Surface_dataset       _surfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_HAND_3_col_hillslope_pftdistribution.nc    Initial_conditions_dataset        finidat_interp_dest.nc     #PFT_physiological_constants_dataset       clm50_params.c240105.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         /./hillslope-soiltest-pfttest.clm2.h0.2011-01.nc    Time_constant_3Dvars      AZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE:PCT_SAND:PCT_CLAY         W   levgrnd                	long_name         coordinate ground levels   units         m         d      AL   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m         P      A�   levlak        	         	long_name         coordinate lake levels     units         m         (      B    levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m               B(   hillslope_distance                 	long_name         hillslope column distance      units         m               B,   hillslope_width                	long_name         hillslope column width     units         m               BD   hillslope_area                 	long_name         hillslope column area      units         m               B\   hillslope_elev                 	long_name         hillslope column elevation     units         m               Bt   hillslope_slope                	long_name         hillslope column slope     units         m               B�   hillslope_aspect               	long_name         hillslope column aspect    units         m               B�   hillslope_index                	long_name         hillslope index             B�   hillslope_cold                 	long_name         hillslope downhill column index             B�   hillslope_colu                 	long_name         hillslope uphill column index               B�   time               	long_name         time   units         days since 2011-01-01 00:00:00     calendar      noleap     bounds        time_bounds             G�   mcdate                 	long_name         current date (YYYYMMDD)             G�   mcsec                  	long_name         current seconds of current date    units         s               G�   mdcur                  	long_name         current day (from base day)             G�   mscur                  	long_name         current seconds of current day              G�   nstep                  	long_name         	time step               G�   time_bounds                   	long_name         history time interval endpoints             G�   date_written                             G�   time_written                             G�   lon                 	long_name         coordinate longitude   units         degrees_east   
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
_FillValue        Gh��y �      `      D    
pfts1d_lat                 	long_name         pft latitude   units         degrees_north      
_FillValue        Gh��y �      `      D`   
pfts1d_ixy                 	long_name         '2d longitude index of corresponding pft    
_FillValue        ����      0      D�   
pfts1d_jxy                 	long_name         &2d latitude index of corresponding pft     
_FillValue        ����      0      D�   	pfts1d_gi                  	long_name         "1d grid index of corresponding pft     
_FillValue        ����      0      E    	pfts1d_li                  	long_name         &1d landunit index of corresponding pft     
_FillValue        ����      0      EP   	pfts1d_ci                  	long_name         $1d column index of corresponding pft   
_FillValue        ����      0      E�   pfts1d_wtgcell                 	long_name         -pft weight relative to corresponding gridcell      
_FillValue        Gh��y �      `      E�   pfts1d_wtlunit                 	long_name         -pft weight relative to corresponding landunit      
_FillValue        Gh��y �      `      F   pfts1d_wtcol               	long_name         +pft weight relative to corresponding column    
_FillValue        Gh��y �      `      Fp   pfts1d_itype_veg               	long_name         pft vegetation type    
_FillValue        ����      0      F�   pfts1d_itype_col               	long_name         'pft column type (see global attributes)    
_FillValue        ����      0      G    pfts1d_itype_lunit                 	long_name         Gpft landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)    
_FillValue        ����      0      G0   pfts1d_active                  	long_name         #true => do computations on this pft    
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                    0      G`   FSAT                  	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            G�   FSNO                  	long_name         "fraction of ground covered by snow     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            G�   H2OSNO                    	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            G�   QICE                  	long_name         ice growth/melt    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            G�   QINFL                     	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            H   QOVER                     	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            H   QRUNOFF                   	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            H    QSNOMELT                  	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            H,   QSOIL                     	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      0      H8   QVEGT                     	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      0      Hh   RAIN                  	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            H�   SNOW                  	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            H�   TREFMNAV                  	long_name         (daily minimum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      0      H�   TREFMXAV                  	long_name         (daily maximum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      0      H�   TSA                   	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      0      I   ZWT                   	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I@   	ZWT_PERCH                     	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            IL   TSOI                     	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��     ,      IX   H2OSOI                       	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      J�   SOILICE                      	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      Kt<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  @      @]�;G�a�@j�7��@d      @�V��%�
@��	� @��     @�     @�T             ?�o��b��@�s���>�MP��>��I�+Ј>ƍH�s�@2ƂH@.���f@�l�:~         ����         ��������C�A�B>@�=xČ?�           @p�1&�x�@G�bM��      @p�1&�x�@G�bM��         ?�            @p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��                                    ?�3��ԣz?ڡK�6�?�h�i�2?�3��ԣz?ڡK�6�?�h�i�2                                    @p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��                                                                                                                                                                                    ?m92��d        ?��j�2]�?�{R�?ʡK�6�                ?ʡK�6�?ñW�U�?ӱW�U�        ?�Au�C?m92��d        ?��j�2]�?�{R�?ʡK�6�                ?ʡK�6�?ñW�U�?ӱW�U�        ?�Au�C?�������        ?�333333?�ffffff?�                      ?�      ?�333334?�333333        ?�������         
            
            
                                                                                                                   D6� 3)5      �      ��@��     @��     04/18/24        16:49:43        >�Q�>�)�>s;�?w��?p��?N�B�9�B�%�A��{@��{@��{@��2N��s��͔8���4Ԅ�4���8�s5#�06�1o6 ��6� 6�+6)�a{@��62[�6�6-��{@��{@��6�74��4$H�{@��4.��4��{@��2�    4%&{@��{@��    4�C2a�{@��V�)4��[4��4�!�7��7��7 C�C�{@��C�G C�E.C�I�{@��{@��C�HEC�G�C�e5{@��C�E$C��{@��C���C��tC�
�{@��{@��C��C��C���{@��C��C��W{@��C���C��vC��{@��{@��C��rC���C��*{@��C��b>�`�>�CR?�@=            C��MC��!C�w�C��}C���C���C��C�zC���C�}�C�r�C���C�ύC��C��C�)C�k�C�i�C���C�TC��]C���C���C�D�C�b�C�:�C���C��kC���C�0�C�ދC��C���C���C�C� }C���C��C��!C���C���C�
1C���C��C��XC��*C��EC��C�C���C�@C�
jC��2C�x#C�FC���C��rC�&�C�y�C���C�8�C�q�C�SsC�N�C�n:C���C�\bC�g�C�r�C�`2C�bC�dC�b�C�b�C�b�?�F�?�BA>8��?l1?l1>B.�?_2S?d*�>9�u?Mj�?Tw�>:��?D�.?>K#>>?Y>`?W��>S?Tz�?Tz�>T��?Tz�?Tz�>D�@?Tz�?Tz�>;��?Tz�?Tz�>,#�?(�k?(�k>>�?(�k?(�k>Q�W>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�wA��%A��&@\B �cB@�ØB�zB��@��A���A��A!?�p?vg�AY$^        A�vW        A�+�        A�81        A}|                                                                                                                                    