CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          '   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 04/18/24 16:56:10   source        #Community Terrestrial Systems Model    hostname      derecho    username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
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
                                                                                                                   E8p 4�      �     )P@��     @�     04/18/24        16:56:10        >�x>��>K�%?{sW?t>r?R�C�<B��A�G�{@��{@��{@��1��۴�譲�5�8���2FU�4ڭ�8�/�45�6R�576O��6��56'2{@��6&Z�6�R6*H{@��{@��6"�:M�2B��{@��3t�q4�n{@��2�:    4��({@��{@��    4H�]���{@��    3�5�3�~�3��6��6��6켖C���{@��C���C��C��-{@��{@��C�̼C��tC��\{@��C�ˆC���{@��C�N�C�NC���{@��{@��C�U�C��}C�]�{@��C�M]C���{@��C���C��DC��k{@��{@��C��iC��@C��{@��C���? ([>�=�@"'o            C�̡C��C�V�C��C�jC�h�C�*zC�>?C��C���C��C��vC��C�T�C��5C�EpC��ZC�<�C��
C��'C���C���C�9�C�XC�6�C��jC��YC���C��3C�D�C��C��C�ܴC�6�C��C�T�C�J�C�6 C��-C�R�C�A�C�HUC�S:C�CC�ʩC�L�C�8�C�L�C�@�C�$C��hC�3�C�	C�(C�(�C��C�o�C�"C��7C���C��C��PC���C�#�C���C�B�C�-�C��GC��eC�>C�t5C��nC�O�C�i*C�z�?�\/?�A�>+��?l1?l1>T�l?ix�?ix�>D\�?e�T?e�T>>q?`Ĝ?A��>N��?D�<?;'�>T*�?8P�?RPq>O�?Ti#?Tz�>T�N?Tz�?Tz�>T�q?Tz�?Tz�>0r ?(�k?(�k>)ږ?(�k?(�k>+wM>��>��>H��>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�w>��>��>T�wA��A�YG@ӐBڨB�j@�6�B:�BB8%t@�|�Bc�jB\��A�=B��BA��B�V?��A�G�?�b�    A��        A��        Bʍ        A���        @�                                                                                                            