CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          '   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 04/16/24 10:38:39   source        #Community Terrestrial Systems Model    hostname      derecho    username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        hillslope-wte-calib-pfts   case_id       hillslope-wte-calib-pfts   Surface_dataset       hsurfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_HAND_3_col_hillslope_wte_pft_calib_fmax0.1006.nc   Initial_conditions_dataset        finidat_interp_dest.nc     #PFT_physiological_constants_dataset       clm50_params.c240105.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         -./hillslope-wte-calib-pfts.clm2.h0.2011-01.nc      Time_constant_3Dvars      AZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE:PCT_SAND:PCT_CLAY           levgrnd                	long_name         coordinate ground levels   units         m         d     �   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m         P        levlak        	         	long_name         coordinate lake levels     units         m         (     l   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m              �   time               	long_name         time   units         days since 2011-01-01 00:00:00     calendar      noleap     bounds        time_bounds            �   mcdate                 	long_name         current date (YYYYMMDD)            �   mcsec                  	long_name         current seconds of current date    units         s              �   mdcur                  	long_name         current day (from base day)            �   mscur                  	long_name         current seconds of current day             �   nstep                  	long_name         	time step              �   time_bounds                   	long_name         history time interval endpoints            �   date_written                            �   time_written                            �   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��           �   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��           �   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��           �   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��           �   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����           �   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����           �   nbedrock                	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����           �   ATM_TOPO                   	long_name         atmospheric surface height     units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   BCDEP                      	long_name         -total BC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                BTRANMN                    	long_name         *daily minimum of transpiration beta factor     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg               DSL                    	long_name         dry surface layer thickness    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               DSTDEP                     	long_name         /total dust deposition (dry+wet) from atmosphere    units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               DSTFLXT                    	long_name         total surface dust emission    units         kg/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               	EFLXBUILD                      	long_name         Cbuilding heat flux from change in interior building air temperature    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               EFLX_DYNBAL                    	long_name         0dynamic land cover change conversion energy flux   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               EFLX_GRND_LAKE                     	long_name         Bnet heat flux into lake/snow surface, excluding light transmission     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               EFLX_LH_TOT                    	long_name         !total latent heat flux [+ to atm]      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                EFLX_LH_TOT_R                      	long_name         Rural total evaporation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   ELAI                   	long_name         !exposed one-sided leaf area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   ERRH2O                     	long_name         total water conservation error     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   	ERRH2OSNO                      	long_name         &imbalance in snow depth (liquid water)     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   ERRSEB                     	long_name         !surface energy conservation error      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   ERRSOI                     	long_name         #soil/lake energy conservation error    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   ERRSOL                     	long_name         "solar radiation conservation error     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   ESAI                   	long_name         !exposed one-sided stem area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   FCEV                   	long_name         canopy evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   FCOV                   	long_name         fractional impermeable area    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            H   FCTR                   	long_name         canopy transpiration   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   FGEV                   	long_name         ground evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   FGR                    	long_name         Oheat flux into soil/snow including snow melt and lake / snow light transmission    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   FGR12                      	long_name         %heat flux between soil layers 1 and 2      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   FH2OSFC                    	long_name         +fraction of ground covered by surface water    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   FIRA                   	long_name         !net infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   FIRA_R                     	long_name         'Rural net infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   FIRE                   	long_name         %emitted infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   FIRE_R                     	long_name         +Rural emitted infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   FLDS                   	long_name         Matmospheric longwave radiation (downscaled for glacier and hillslope columns)      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   FPSN                   	long_name         photosynthesis     units         umol m-2 s-1   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   FSA                    	long_name         absorbed solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   FSAT                   	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            |   FSDS                   	long_name         Satmospheric incident solar radiation (downscaled for glacier and hillslope columns)    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSND                     	long_name         #direct nir incident solar radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSNDLN                   	long_name         1direct nir incident solar radiation at local noon      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSNI                     	long_name         $diffuse nir incident solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSVD                     	long_name         #direct vis incident solar radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSVDLN                   	long_name         1direct vis incident solar radiation at local noon      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSVI                     	long_name         $diffuse vis incident solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSVILN                   	long_name         2diffuse vis incident solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDS_from_atm                      	long_name         Oatmospheric incident solar radiation received from atmosphere (pre-downscaling)    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSH                    	long_name         Ssensible heat not including correction for land use change and rain/snow conversion    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSH_G                      	long_name         sensible heat from ground      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSH_PRECIP_CONVERSION                      	long_name         ;Sensible heat flux from conversion of rain/snow atm forcing    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSH_R                      	long_name         Rural sensible heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSH_RUNOFF_ICE_TO_LIQ                      	long_name         Dsensible heat flux generated from conversion of ice runoff to liquid   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSH_TO_COUPLER                     	long_name         �sensible heat sent to coupler (includes corrections for land use change, rain/snow conversion and conversion of ice runoff to liquid)      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSH_V                      	long_name         sensible heat from veg     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSM                    	long_name         snow melt heat flux    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSNO                   	long_name         "fraction of ground covered by snow     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSNO_EFF                   	long_name         ,effective fraction of ground covered by snow   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSR                    	long_name         reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSRND                      	long_name         $direct nir reflected solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSRNDLN                    	long_name         2direct nir reflected solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSRNI                      	long_name         %diffuse nir reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSRVD                      	long_name         $direct vis reflected solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSRVDLN                    	long_name         2direct vis reflected solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSRVI                      	long_name         %diffuse vis reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   GSSHA                      	long_name          shaded leaf stomatal conductance   units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   GSSHALN                    	long_name         .shaded leaf stomatal conductance at local noon     units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   GSSUN                      	long_name          sunlit leaf stomatal conductance   units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   GSSUNLN                    	long_name         .sunlit leaf stomatal conductance at local noon     units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   H2OCAN                     	long_name         intercepted water      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   H2OSFC                     	long_name         surface water depth    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   H2OSNO                     	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                
H2OSNO_TOP                     	long_name         mass of snow in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HEAT_CONTENT1                      	long_name         #initial gridcell total heat content    units         J/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HEAT_FROM_AC                   	long_name         Lsensible heat flux put into canyon due to heat removed from air conditioning   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HIA                    	long_name         2 m NWS Heat Index     units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HIA_R                      	long_name         Rural 2 m NWS Heat Index   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HIA_U                      	long_name         Urban 2 m NWS Heat Index   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HUMIDEX                    	long_name         2 m Humidex    units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               	HUMIDEX_R                      	long_name         Rural 2 m Humidex      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                	HUMIDEX_U                      	long_name         Urban 2 m Humidex      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   ICE_CONTENT1                   	long_name         "initial gridcell total ice content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   IWUELN                     	long_name         )local noon intrinsic water use efficiency      units         umolCO2/molH2O     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   JMX25T                     	long_name         canopy profile of jmax     units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   Jmx25Z                     	long_name         Bmaximum rate of electron transport at 25 Celcius for canopy layers     units         umol electrons/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   LAISHA                     	long_name          shaded projected leaf area index   units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   LAISUN                     	long_name          sunlit projected leaf area index   units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   LAKEICEFRAC_SURF                   	long_name         $surface lake layer ice mass fraction   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   LAKEICETHICK                   	long_name         @thickness of lake ice (including physical expansion on freezing)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   LIQCAN                     	long_name         intercepted liquid water   units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   LIQUID_CONTENT1                    	long_name         "initial gridcell total liq content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   LNC                    	long_name         leaf N concentration   units         gN leaf/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   MEG_acetaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   MEG_acetic_acid                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   MEG_acetone                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   MEG_carene_3                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   MEG_ethanol                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   MEG_formaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   MEG_isoprene                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   MEG_methanol                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   MEG_pinene_a                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   MEG_thujene_a                      	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   OCDEP                      	long_name         -total OC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   PARVEGLN                   	long_name         (absorbed par by vegetation at local noon   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   PBOT                   	long_name         Natmospheric pressure at surface (downscaled for glacier and hillslope columns)     units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   PCO2                   	long_name         #atmospheric partial pressure of CO2    units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   Q2M                    	long_name         2m specific humidity   units         kg/kg      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QBOT                   	long_name         Hatmospheric specific humidity (downscaled to columns in glacier regions)   units         kg/kg      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QDRAI                      	long_name         sub-surface drainage   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QDRAI_PERCH                    	long_name         perched wt drainage    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QDRAI_XS                   	long_name         saturation excess drainage     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLOOD                     	long_name         runoff from river flooding     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_EVAP_TOT                      	long_name         -qflx_evap_soi + qflx_evap_can + qflx_tran_veg      units         
kg m-2 s-1     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_ICE_DYNBAL                    	long_name         4ice dynamic land cover change conversion runoff flux   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_LIQDEW_TO_TOP_LAYER                   	long_name         >rate of liquid water deposited on top soil or snow layer (dew)     units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_LIQEVAP_FROM_TOP_LAYER                    	long_name         ;rate of liquid water evaporated from top soil or snow layer    units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_LIQ_DYNBAL                    	long_name         4liq dynamic land cover change conversion runoff flux   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_SNOW_DRAIN                    	long_name         drainage from snow pack    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_SNOW_DRAIN_ICE                    	long_name         1drainage from snow pack melt (ice landunits only)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            �   QFLX_SOLIDDEW_TO_TOP_LAYER                     	long_name         ?rate of solid water deposited on top soil or snow layer (frost)    units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_SOLIDEVAP_FROM_TOP_LAYER                      	long_name         zrate of ice evaporated from top soil or snow layer (sublimation) (also includes bare ice sublimation from glacier columns)     units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QH2OSFC                    	long_name         surface water runoff   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QHR                    	long_name         hydraulic redistribution   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   QICE                   	long_name         ice growth/melt    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            �   QICE_FRZ                   	long_name         
ice growth     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            �   	QICE_MELT                      	long_name         ice melt   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            �   QINFL                      	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QINTR                      	long_name         interception   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QIRRIG_FROM_GW_CONFINED                    	long_name         3water added through confined groundwater irrigation    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QIRRIG_FROM_GW_UNCONFINED                      	long_name         5water added through unconfined groundwater irrigation      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QIRRIG_FROM_SURFACE                    	long_name         ,water added through surface water irrigation   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QLATFLOWOUT                    	long_name         hillcol lateral outflow    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg             �   QOVER                      	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QRGWL                      	long_name         isurface runoff at glaciers (liquid only), wetlands, lakes; also includes melted ice runoff from QSNWCPICE      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QRUNOFF                    	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QRUNOFF_ICE                    	long_name         Btotal liquid runoff not incl corret for LULCC (ice landunits only)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice                QRUNOFF_ICE_TO_COUPLER                     	long_name         Ktotal ice runoff sent to coupler (includes corrections for land use change)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QRUNOFF_TO_COUPLER                     	long_name         Ntotal liquid runoff sent to coupler (includes corrections for land use change)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               	QSNOCPLIQ                      	long_name         Rexcess liquid h2o due to snow capping not including correction for land use change     units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QSNOEVAP                   	long_name         Nevaporation from snow (only when snl<0, otherwise it is equal to qflx_ev_soil)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QSNOFRZ                    	long_name         $column-integrated snow freezing rate   units         kg/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QSNOFRZ_ICE                    	long_name         9column-integrated snow freezing rate (ice landunits only)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice               QSNOMELT                   	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QSNOMELT_ICE                   	long_name         snow melt (ice landunits only)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice                
QSNOUNLOAD                     	long_name         canopy snow unloading      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   QSNO_TEMPUNLOAD                    	long_name         canopy snow temp unloading     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   QSNO_WINDUNLOAD                    	long_name         canopy snow wind unloading     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   	QSNWCPICE                      	long_name         Qexcess solid h2o due to snow capping not including correction for land use change      units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   QSOIL                      	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   	QSOIL_ICE                      	long_name         'Ground evaporation (ice landunits only)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            8   QVEGE                      	long_name         canopy evaporation     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   QVEGT                      	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   RAIN                   	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   RAIN_FROM_ATM                      	long_name         >atmospheric rain received from atmosphere (pre-repartitioning)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   RH2M                   	long_name         2m relative humidity   units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   RSSHA                      	long_name         shaded leaf stomatal resistance    units         s/m    cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            P   RSSUN                      	long_name         sunlit leaf stomatal resistance    units         s/m    cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            T   SABG                   	long_name         solar rad absorbed by ground   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   SABG_PEN                   	long_name         2Rural solar rad penetrating top soil or snow layer     units         watt/m^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   SABV                   	long_name         solar rad absorbed by veg      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   SNOBCMCL                   	long_name         mass of BC in snow column      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   SNOBCMSL                   	long_name         mass of BC in top snow layer   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   SNOCAN                     	long_name         intercepted snow   units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   	SNODSTMCL                      	long_name         mass of dust in snow column    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   	SNODSTMSL                      	long_name         mass of dust in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   SNOFSRND                   	long_name         .direct nir reflected solar radiation from snow     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   SNOFSRNI                   	long_name         /diffuse nir reflected solar radiation from snow    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   SNOFSRVD                   	long_name         .direct vis reflected solar radiation from snow     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOFSRVI                   	long_name         /diffuse vis reflected solar radiation from snow    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	SNOINTABS                      	long_name         8Fraction of incoming solar absorbed by lower snow layers   units         -      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOMELT_ACCUM                      	long_name         accumulated snow melt for z0   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOOCMCL                   	long_name         mass of OC in snow column      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOOCMSL                   	long_name         mass of OC in top snow layer   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	SNOTXMASS                      	long_name         ksnow temperature times layer mass, layer sum; to get mass-weighted temperature, divide by (SNOWICE+SNOWLIQ)    units         K kg/m2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOW                   	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOWDP                     	long_name         gridcell mean snow height      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOWICE                    	long_name         snow ice   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOWLIQ                    	long_name         snow liquid water      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   
SNOW_DEPTH                     	long_name          snow height of snow covered area   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOW_FROM_ATM                      	long_name         >atmospheric snow received from atmosphere (pre-repartitioning)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOW_PERSISTENCE                   	long_name         BLength of time of continuous snow cover (nat. veg. landunits only)     units         seconds    cell_methods      time: instantaneous    
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg             �   
SNOW_SINKS                     	long_name         snow sinks (liquid water)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOW_SOURCES                   	long_name         snow sources (liquid water)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	SOILRESIS                      	long_name         soil resistance to evaporation     units         s/m    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SOILWATER_10CM                     	long_name         @soil liquid water + ice in top 10cm of soil (veg landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   SWBGT                      	long_name         !2 m Simplified Wetbulb Globe Temp      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SWBGT_R                    	long_name         'Rural 2 m Simplified Wetbulb Globe Temp    units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SWBGT_U                    	long_name         'Urban 2 m Simplified Wetbulb Globe Temp    units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TAUX                   	long_name         zonal surface stress   units         kg/m/s^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TAUY                   	long_name         meridional surface stress      units         kg/m/s^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TBOT                   	long_name         Jatmospheric air temperature (downscaled for glacier and hillslope columns)     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TBUILD                     	long_name         'internal urban building air temperature    units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TG                     	long_name         ground temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TH2OSFC                    	long_name         surface water temperature      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   THBOT                      	long_name         Tatmospheric air potential temperature (downscaled for glacier and hillslope columns)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TKE1                   	long_name         (top lake level eddy thermal conductivity   units         W/(mK)     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TLAI                   	long_name         total projected leaf area index    units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   
TOTSOILICE                     	long_name         /vertically summed soil ice (veg landunits only)    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   
TOTSOILLIQ                     	long_name         8vertically summed soil liquid water (veg landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   TPU25T                     	long_name         canopy profile of tpu      units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                TREFMNAV                   	long_name         (daily minimum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TREFMXAV                   	long_name         (daily maximum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TSA                    	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TSAI                   	long_name         total projected stem area index    units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TSKIN                      	long_name         skin temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TSL                    	long_name         Rtemperature of near-surface soil layer (natural vegetated and crop landunits only)     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg               	TSOI_10CM                      	long_name         $soil temperature in top 10cm of soil   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TV                     	long_name         vegetation temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                TWS                    	long_name         total water storage    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   U10                    	long_name         	10-m wind      units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   U10_DUST                   	long_name         10-m wind for dust model   units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   URBAN_AC                   	long_name         urban air conditioning flux    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   
URBAN_HEAT                     	long_name         urban heating flux     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   VCMX25T                    	long_name         canopy profile of vcmax25      units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   VENTILATION                    	long_name         ,sensible heat flux from building ventilation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   VOLR                   	long_name         !river channel total water storage      units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   VOLRMCH                    	long_name         (river channel main channel water storage   units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   VOLUMETRIC_STREAMFLOW                      	long_name         $volumetric streamflow from hillslope   units         m3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg             H   VPD2M                      	long_name         2m vapor pressure deficit      units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   VPD_CAN                    	long_name         canopy vapor pressure deficit      units         kPa    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   Vcmx25Z                    	long_name         1canopy profile of vcmax25 predicted by LUNA model      units         	umol/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   	WASTEHEAT                      	long_name         Csensible heat flux from heating/cooling sources of urban waste heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   WBT                    	long_name         2 m Stull Wet Bulb     units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   WBT_R                      	long_name         Rural 2 m Stull Wet Bulb   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   WBT_U                      	long_name         Urban 2 m Stull Wet Bulb   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   WIND                   	long_name         #atmospheric wind velocity magnitude    units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   ZBOT                   	long_name         atmospheric reference height   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   ZWT                    	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            p   	ZWT_PERCH                      	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            t   PCT_CFT                       	long_name         #% of each crop on the crop landunit    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   PCT_GLC_MEC                       	long_name         5% of each GLC elevation class on the glacier landunit      units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       (     �   SMP                       	long_name         Asoil matric potential (natural vegetated and crop landunits only)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d     �   TSOI                      	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d        TSOI_ICE                      	long_name         %soil temperature (ice landunits only)      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice       d     p   TLAKE            	             	long_name         lake temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       (     �   H2OSOI                        	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     �   SOILICE                       	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     L   SOILLIQ                       	long_name         =soil liquid water (natural vegetated and crop landunits only)      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     �   PCT_LANDUNIT                      	long_name         % of each landunit on grid cell    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       $     �   PCT_NAT_PFT                       	long_name         =% of each PFT on the natural vegetation (i.e., soil) landunit      units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       <        VEGWP                         	long_name         Fvegetation water matric potential for sun/sha canopy,xyl,root segments     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   VEGWPLN                       	long_name         Kvegetation water matric potential for sun/sha canopy,xyl,root at local noon    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   VEGWPPD                       	long_name         Epredawn vegetation water matric potential for sun/sha canopy,xyl,root      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  C�A�B>@�=xČ?�           D�  3�e      !     V0@�     @��     04/16/24        10:38:39        C�Ws+nV?�v<�D/+�u            {@��@��@��;6��(���#�@���S�-�T
��'>�"Ҽ�w�=�o;�3�@�VC����t�i<���@�Z5@�Z5C���C���C��f:���@�X=�oAc��@c?�A���@d`
@��AZu�@��|A�w8Ac��@C^?d{@�!�N@C^    @C�@	�@b�m?t!	?t!	@�x?��@�`P?��L?�G�A	��@Hi�G)�G��GM՞G���>�NSB���C19�@i��N��y    ��LT��LT{@��� ��� ��{@��C��2A�݁BcBK�c;3�:C��{@��{@��=�&E���?��t%�8#��v'	È&y6%�8#���&�B'e�&r��$/��,ô_A#��G�<BXy;0��;/�\�f;�            5�%&    4��2��    6��{@��2b�75�fS6|�ôDv�{@��{@��{@�γ��X6��*            56���    6u�^{@��    6u�^    5��87~�{@��72�{@��64�6�s5D�_    5��{@�α�z 0�46�A+6�=3B���        @�3?�FU@�6�5���45	�>�܏6c�4�rx@	�@mc@�@��	=��5;���7^[�5��vG<6�7c`t>뵬C-'2@�SF>�W�7cbpJ_�F7�C7�١A�1�B�ve@6�y@6�y{@�λ1�d�1�dC�(B{@��C�A2C���C�(B{@��=�"B��RE��@!�"C��4C���C�.?9K�C�QBC�jC�w�C�inE���?�=�?uӉ        Aq��            9iH�A��=|�A�(�    ���m���m{@��?k`A   >���                                                    �!���!��<�U��%r�_�G���!�����)Fô�düQ��2���2���2���2���2���2���2���2��̾� ̾� ̾� ̾� ̾� C�jC�r>C��*C��UC�${C���C�SC��C��C��RC�+�C�7
C�3�C�*�C�}C��C���C��HC��C��=C�C�$@C�G�C�b�C�q{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��?���?l1?e��?Ϻ?X9$?Y��?Tz�?Tz�?Tz�?Tz�?(�k?(�k>��>��>��>��>��>��>��>��A�ۍA��yBE�?�                                                                ?s�@"lAJQ�B;��Bʵ�C  C&  CG33ChffC���Cl�fC���CVC�YC��qC��cC�RTC�FC��8D
��B�                                                              B&p        @��        BMh    �W��U$�STç��������yù���/��/��/�ã�B