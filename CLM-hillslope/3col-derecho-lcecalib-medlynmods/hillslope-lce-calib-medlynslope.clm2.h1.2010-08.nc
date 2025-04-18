CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          '   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 03/31/25 16:54:07   source        #Community Terrestrial Systems Model    hostname      derecho    username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        hillslope-lce-calib-medlynslope    case_id       hillslope-lce-calib-medlynslope    Surface_dataset       Zsurfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_HAND_3_col_hillslope_medlynmods.nc     Initial_conditions_dataset        finidat_interp_dest.nc     #PFT_physiological_constants_dataset       ,clm50_params.c240105_hillslope_medlynmods.nc   ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         4./hillslope-lce-calib-medlynslope.clm2.h0.2004-01.nc   Time_constant_3Dvars      AZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE:PCT_SAND:PCT_CLAY         X   levgrnd                	long_name         coordinate ground levels   units         m         d      B@   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m         P      B�   levlak        	         	long_name         coordinate lake levels     units         m         (      B�   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m               C   hillslope_distance                 	long_name         hillslope column distance      units         m               C    hillslope_width                	long_name         hillslope column width     units         m               C8   hillslope_area                 	long_name         hillslope column area      units         m               CP   hillslope_elev                 	long_name         hillslope column elevation     units         m               Ch   hillslope_slope                	long_name         hillslope column slope     units         m               C�   hillslope_aspect               	long_name         hillslope column aspect    units         m               C�   hillslope_index                	long_name         hillslope index             C�   hillslope_cold                 	long_name         hillslope downhill column index             C�   hillslope_colu                 	long_name         hillslope uphill column index               C�   time               	long_name         time   units         days since 2004-01-01 00:00:00     calendar      noleap     bounds        time_bounds             Ih   mcdate                 	long_name         current date (YYYYMMDD)             Il   mcsec                  	long_name         current seconds of current date    units         s               Ip   mdcur                  	long_name         current day (from base day)             It   mscur                  	long_name         current seconds of current day              Ix   nstep                  	long_name         	time step               I|   time_bounds                   	long_name         history time interval endpoints             I�   date_written                             I�   time_written                             I�   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��            C�   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��            C�   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��            C�   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��            C�   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����            C�   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����            C�   nbedrock                	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����            C�   
grid1d_lon                 	long_name         gridcell longitude     units         degrees_east   
_FillValue        Gh��y �            C�   
grid1d_lat                 	long_name         gridcell latitude      units         degrees_north      
_FillValue        Gh��y �            C�   
grid1d_ixy                 	long_name         ,2d longitude index of corresponding gridcell   
_FillValue        ����            D    
grid1d_jxy                 	long_name         +2d latitude index of corresponding gridcell    
_FillValue        ����            D   
land1d_lon                 	long_name         landunit longitude     units         degrees_east   
_FillValue        Gh��y �            D   
land1d_lat                 	long_name         landunit latitude      units         degrees_north      
_FillValue        Gh��y �            D   
land1d_ixy                 	long_name         ,2d longitude index of corresponding landunit   
_FillValue        ����            D   
land1d_jxy                 	long_name         +2d latitude index of corresponding landunit    
_FillValue        ����            D   	land1d_gi                  	long_name         '1d grid index of corresponding landunit    
_FillValue        ����            D    land1d_wtgcell                 	long_name         2landunit weight relative to corresponding gridcell     
_FillValue        Gh��y �            D$   land1d_ityplunit               	long_name         Clandunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)    
_FillValue        ����            D,   land1d_active                  	long_name         (true => do computations on this landunit   
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          D0   
cols1d_lon                 	long_name         column longitude   units         degrees_east   
_FillValue        Gh��y �            D4   
cols1d_lat                 	long_name         column latitude    units         degrees_north      
_FillValue        Gh��y �            DL   
cols1d_ixy                 	long_name         *2d longitude index of corresponding column     
_FillValue        ����            Dd   
cols1d_jxy                 	long_name         )2d latitude index of corresponding column      
_FillValue        ����            Dp   	cols1d_gi                  	long_name         %1d grid index of corresponding column      
_FillValue        ����            D|   	cols1d_li                  	long_name         )1d landunit index of corresponding column      
_FillValue        ����            D�   cols1d_wtgcell                 	long_name         0column weight relative to corresponding gridcell   
_FillValue        Gh��y �            D�   cols1d_wtlunit                 	long_name         0column weight relative to corresponding landunit   
_FillValue        Gh��y �            D�   cols1d_itype_col               	long_name         #column type (see global attributes)    
_FillValue        ����            D�   cols1d_itype_lunit                 	long_name         Jcolumn landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)     
_FillValue        ����            D�   cols1d_active                  	long_name         &true => do computations on this column     
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                          D�   cols1d_nbedrock                	long_name         column bedrock depth index     
_FillValue        ����            D�   
pfts1d_lon                 	long_name         pft longitude      units         degrees_east   
_FillValue        Gh��y �      x      D�   
pfts1d_lat                 	long_name         pft latitude   units         degrees_north      
_FillValue        Gh��y �      x      El   
pfts1d_ixy                 	long_name         '2d longitude index of corresponding pft    
_FillValue        ����      <      E�   
pfts1d_jxy                 	long_name         &2d latitude index of corresponding pft     
_FillValue        ����      <      F    	pfts1d_gi                  	long_name         "1d grid index of corresponding pft     
_FillValue        ����      <      F\   	pfts1d_li                  	long_name         &1d landunit index of corresponding pft     
_FillValue        ����      <      F�   	pfts1d_ci                  	long_name         $1d column index of corresponding pft   
_FillValue        ����      <      F�   pfts1d_wtgcell                 	long_name         -pft weight relative to corresponding gridcell      
_FillValue        Gh��y �      x      G   pfts1d_wtlunit                 	long_name         -pft weight relative to corresponding landunit      
_FillValue        Gh��y �      x      G�   pfts1d_wtcol               	long_name         +pft weight relative to corresponding column    
_FillValue        Gh��y �      x      H    pfts1d_itype_veg               	long_name         pft vegetation type    
_FillValue        ����      <      Hx   pfts1d_itype_col               	long_name         'pft column type (see global attributes)    
_FillValue        ����      <      H�   pfts1d_itype_lunit                 	long_name         Gpft landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)    
_FillValue        ����      <      H�   pfts1d_active                  	long_name         #true => do computations on this pft    
_FillValue               flag_values                 flag_meanings         
FALSE TRUE     valid_range                    <      I,   FSAT                  	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I�   FSNO                  	long_name         "fraction of ground covered by snow     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I�   H2OSNO                    	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I�   QICE                  	long_name         ice growth/melt    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I�   QINFL                     	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I�   QOVER                     	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I�   QRUNOFF                   	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            I�   QSNOMELT                  	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            J   QSOIL                     	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      <      J   QVEGT                     	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      <      JL   RAIN                  	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            J�   SNOW                  	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            J�   TREFMNAV                  	long_name         (daily minimum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      <      J�   TREFMXAV                  	long_name         (daily maximum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      <      J�   TSA                   	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      <      K   ZWT                   	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            KT   	ZWT_PERCH                     	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��            K`   TSOI                     	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��     ,      Kl   watsat                       	long_name         water saturated    units         m^3/m^3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��     ,      L�   H2OSOI                       	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      M�   SOILICE                      	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��      �      N�<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  @      @]�;G�a�@j�7��@d      @�V��%�
@��	� @��     @�     @�T             ?�o��b��@�s���>�MP��>��I�+Ј>ƍH�s�@2ƂH@.���f@�l�:~         ����         ��������C�A�B>@�=xČ?�           @p�1&�x�@G�bM��      @p�1&�x�@G�bM��         ?�            @p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��                                    ?�3��ԣz?ڡK�6�?�h�i�2?�3��ԣz?ڡK�6�?�h�i�2                                    @p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@p�1&�x�@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��@G�bM��                                                                                                                                                                                                                                 ?�92��]?�92��[        ?�a�,��T?��j�2]�?¤z@Q?���cۉ�        ?�M�B�[�?���cۉ�?�Au�C        ?� �$?��G��x?�Au�C?�92��]?�92��[        ?�a�,��T?��j�2]�?¤z@Q?���cۉ�        ?�M�B�[�?���cۉ�?�Au�C        ?� �$?��G��x?�Au�C?ə�����?ə�����        ?�������?�333333?�ffffff?�333332        ?ə�����?�333333?�������        ?ᙙ����?�fffffe?�������                                                                                                                                                                                       E 2�%      	�     �0@��     @�     03/31/25        16:54:07                                            {@��{@��{@��=\z%7��l7��#            �*͸5 �6+��            7��98d�{@��7��$6oN7�67>�{@��7�69��5��!{@��5�kP5��54�7�V6��C{@��6�rd8	4�8&�7��{@��7�{8�<8�{@��8��7��8�17�Ȍ7��7�/�            C�neC�v�{@��C�g�C�L�C�[�C�e�{@��C�ZPC�B�C��{@��C�ahC�x�C�R�C�6C�L�{@��C� xC�.�C�ZMC�z�{@��C�L�C�5�C�Ck{@��C�D-C�1�C�.�C��+C���{@��C���C��%C��yC���{@��C���C��PC���{@��C��SC��1C�� <�S�>��?�V@�R@�R@�RC��C���C�C�$C���C��C�2C�pC��C��C�B�C���C���C�C���C�7�C���C��FC���C���C���C�¹C�4�C���C��8C�IBC�C�C��yC�T�C���C�?�C���C��\C��C�B�C�7�C���C���C��C��DC��tC�LkC�j�C��C���C�\�C���C�*C�Y�C���C���C�]�C�}9C�C�d�C�zC���C�jC�w�C�ogC�l�C�tC�*+C�j�C�k�C�C�hjC�d�C��C�m4C�j�C���C�pC�o�C���?m�h?m�h>[2 ?l1?l1>[2 ?ix�?ix�>[2 ?e�T?e�T>[2 ?`Ĝ?`Ĝ>[2 ?Y��?Y��>[2 ?Tz�?Tz�>[2 ?Tz�?Tz�>[2 ?Tz�?Tz�>[2 ?Tz�?Tz�>[2 ?(�k?(�k>[2 ?(�k?(�k>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>y�>��>��>y�>��>��>y�>��>��>y�>��>��>y�?m�h>��->V*?l1>��N>`+?ix�>�U#>�5?e�T?l�>8&?`Ĝ?K�G> h�?Y��?Y��>"t�?Tz�?Tz�>%A�?Tz�?Tz�>(��?Tz�?Tz�>,А?Tz�?Tz�>2bV?(�k?(�k>;Rn?(�k?(�k>Pe�>��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2 >��>��>[2                                                                                                                                                                                                                                                 