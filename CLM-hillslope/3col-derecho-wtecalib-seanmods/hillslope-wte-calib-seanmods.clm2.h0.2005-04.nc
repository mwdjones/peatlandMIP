CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          '   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 08/20/24 13:53:04   source        #Community Terrestrial Systems Model    hostname      derecho    username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        hillslope-wte-calib-seanmods   case_id       hillslope-wte-calib-seanmods   Surface_dataset       Xsurfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_HAND_3_col_hillslope_seanmods.nc   Initial_conditions_dataset        finidat_interp_dest.nc     #PFT_physiological_constants_dataset       *clm50_params.c240105_hillslope_seanmods.nc     ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         1./hillslope-wte-calib-seanmods.clm2.h0.2004-01.nc      Time_constant_3Dvars      AZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE:PCT_SAND:PCT_CLAY           levgrnd                	long_name         coordinate ground levels   units         m         d     �   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m         P     ,   levlak        	         	long_name         coordinate lake levels     units         m         (     |   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m              �   time               	long_name         time   units         days since 2004-01-01 00:00:00     calendar      noleap     bounds        time_bounds            �   mcdate                 	long_name         current date (YYYYMMDD)            �   mcsec                  	long_name         current seconds of current date    units         s              �   mdcur                  	long_name         current day (from base day)            �   mscur                  	long_name         current seconds of current day             �   nstep                  	long_name         	time step              �   time_bounds                   	long_name         history time interval endpoints            �   date_written                            �   time_written                            �   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��           �   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��           �   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��           �   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��           �   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����           �   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����           �   nbedrock                	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����           �   ATM_TOPO                   	long_name         atmospheric surface height     units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               BCDEP                      	long_name         -total BC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               BTRANMN                    	long_name         *daily minimum of transpiration beta factor     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg               DSL                    	long_name         dry surface layer thickness    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               DSTDEP                     	long_name         /total dust deposition (dry+wet) from atmosphere    units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               DSTFLXT                    	long_name         total surface dust emission    units         kg/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                	EFLXBUILD                      	long_name         Cbuilding heat flux from change in interior building air temperature    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   EFLX_DYNBAL                    	long_name         0dynamic land cover change conversion energy flux   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   EFLX_GRND_LAKE                     	long_name         Bnet heat flux into lake/snow surface, excluding light transmission     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   EFLX_LH_TOT                    	long_name         !total latent heat flux [+ to atm]      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   EFLX_LH_TOT_R                      	long_name         Rural total evaporation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   ELAI                   	long_name         !exposed one-sided leaf area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   ERRH2O                     	long_name         total water conservation error     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   	ERRH2OSNO                      	long_name         &imbalance in snow depth (liquid water)     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   ERRSEB                     	long_name         !surface energy conservation error      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   ERRSOI                     	long_name         #soil/lake energy conservation error    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   ERRSOL                     	long_name         "solar radiation conservation error     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   ESAI                   	long_name         !exposed one-sided stem area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   FCEV                   	long_name         canopy evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   FCOV                   	long_name         fractional impermeable area    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            X   FCTR                   	long_name         canopy transpiration   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   FGEV                   	long_name         ground evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   FGR                    	long_name         Oheat flux into soil/snow including snow melt and lake / snow light transmission    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   FGR12                      	long_name         %heat flux between soil layers 1 and 2      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   FH2OSFC                    	long_name         +fraction of ground covered by surface water    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   FIRA                   	long_name         !net infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   FIRA_R                     	long_name         'Rural net infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   FIRE                   	long_name         %emitted infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   FIRE_R                     	long_name         +Rural emitted infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   FLDS                   	long_name         Matmospheric longwave radiation (downscaled for glacier and hillslope columns)      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FPSN                   	long_name         photosynthesis     units         umol m-2 s-1   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSA                    	long_name         absorbed solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSAT                   	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   FSDS                   	long_name         Satmospheric incident solar radiation (downscaled for glacier and hillslope columns)    units         W/m^2      cell_methods      
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
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                GSSUNLN                    	long_name         .sunlit leaf stomatal conductance at local noon     units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               H2OCAN                     	long_name         intercepted water      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               H2OSFC                     	long_name         surface water depth    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               H2OSNO                     	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               
H2OSNO_TOP                     	long_name         mass of snow in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HEAT_CONTENT1                      	long_name         #initial gridcell total heat content    units         J/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HEAT_FROM_AC                   	long_name         Lsensible heat flux put into canyon due to heat removed from air conditioning   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               HIA                    	long_name         2 m NWS Heat Index     units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                HIA_R                      	long_name         Rural 2 m NWS Heat Index   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   HIA_U                      	long_name         Urban 2 m NWS Heat Index   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   HUMIDEX                    	long_name         2 m Humidex    units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   	HUMIDEX_R                      	long_name         Rural 2 m Humidex      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   	HUMIDEX_U                      	long_name         Urban 2 m Humidex      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   ICE_CONTENT1                   	long_name         "initial gridcell total ice content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   IWUELN                     	long_name         )local noon intrinsic water use efficiency      units         umolCO2/molH2O     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   JMX25T                     	long_name         canopy profile of jmax     units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   Jmx25Z                     	long_name         Bmaximum rate of electron transport at 25 Celcius for canopy layers     units         umol electrons/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   LAISHA                     	long_name          shaded projected leaf area index   units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   LAISUN                     	long_name          sunlit projected leaf area index   units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   LAKEICEFRAC_SURF                   	long_name         $surface lake layer ice mass fraction   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   LAKEICETHICK                   	long_name         @thickness of lake ice (including physical expansion on freezing)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   LIQCAN                     	long_name         intercepted liquid water   units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   LIQUID_CONTENT1                    	long_name         "initial gridcell total liq content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   LNC                    	long_name         leaf N concentration   units         gN leaf/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   MEG_acetaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   MEG_acetic_acid                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   MEG_acetone                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   MEG_carene_3                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   MEG_ethanol                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   MEG_formaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   MEG_isoprene                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   MEG_methanol                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   MEG_pinene_a                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   MEG_thujene_a                      	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   OCDEP                      	long_name         -total OC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   PARVEGLN                   	long_name         (absorbed par by vegetation at local noon   units         W/m^2      cell_methods      
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
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg                 QOVER                      	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QRGWL                      	long_name         isurface runoff at glaciers (liquid only), wetlands, lakes; also includes melted ice runoff from QSNWCPICE      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QRUNOFF                    	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QRUNOFF_ICE                    	long_name         Btotal liquid runoff not incl corret for LULCC (ice landunits only)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice               QRUNOFF_ICE_TO_COUPLER                     	long_name         Ktotal ice runoff sent to coupler (includes corrections for land use change)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QRUNOFF_TO_COUPLER                     	long_name         Ntotal liquid runoff sent to coupler (includes corrections for land use change)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               	QSNOCPLIQ                      	long_name         Rexcess liquid h2o due to snow capping not including correction for land use change     units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QSNOEVAP                   	long_name         Nevaporation from snow (only when snl<0, otherwise it is equal to qflx_ev_soil)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                QSNOFRZ                    	long_name         $column-integrated snow freezing rate   units         kg/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   QSNOFRZ_ICE                    	long_name         9column-integrated snow freezing rate (ice landunits only)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            (   QSNOMELT                   	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   QSNOMELT_ICE                   	long_name         snow melt (ice landunits only)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            0   
QSNOUNLOAD                     	long_name         canopy snow unloading      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   QSNO_TEMPUNLOAD                    	long_name         canopy snow temp unloading     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   QSNO_WINDUNLOAD                    	long_name         canopy snow wind unloading     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   	QSNWCPICE                      	long_name         Qexcess solid h2o due to snow capping not including correction for land use change      units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   QSOIL                      	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   	QSOIL_ICE                      	long_name         'Ground evaporation (ice landunits only)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            H   QVEGE                      	long_name         canopy evaporation     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   QVEGT                      	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   RAIN                   	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   RAIN_FROM_ATM                      	long_name         >atmospheric rain received from atmosphere (pre-repartitioning)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   RH2M                   	long_name         2m relative humidity   units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   RSSHA                      	long_name         shaded leaf stomatal resistance    units         s/m    cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            `   RSSUN                      	long_name         sunlit leaf stomatal resistance    units         s/m    cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            d   SABG                   	long_name         solar rad absorbed by ground   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   SABG_PEN                   	long_name         2Rural solar rad penetrating top soil or snow layer     units         watt/m^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   SABV                   	long_name         solar rad absorbed by veg      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   SNOBCMCL                   	long_name         mass of BC in snow column      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   SNOBCMSL                   	long_name         mass of BC in top snow layer   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   SNOCAN                     	long_name         intercepted snow   units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   	SNODSTMCL                      	long_name         mass of dust in snow column    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	SNODSTMSL                      	long_name         mass of dust in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOFSRND                   	long_name         .direct nir reflected solar radiation from snow     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOFSRNI                   	long_name         /diffuse nir reflected solar radiation from snow    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SNOFSRVD                   	long_name         .direct vis reflected solar radiation from snow     units         W/m^2      cell_methods      
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
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                TLAI                   	long_name         total projected leaf area index    units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               
TOTSOILICE                     	long_name         /vertically summed soil ice (veg landunits only)    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg               
TOTSOILLIQ                     	long_name         8vertically summed soil liquid water (veg landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg               TPU25T                     	long_name         canopy profile of tpu      units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TREFMNAV                   	long_name         (daily minimum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TREFMXAV                   	long_name         (daily maximum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TSA                    	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               TSAI                   	long_name         total projected stem area index    units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                TSKIN                      	long_name         skin temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   TSL                    	long_name         Rtemperature of near-surface soil layer (natural vegetated and crop landunits only)     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            (   	TSOI_10CM                      	long_name         $soil temperature in top 10cm of soil   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   TV                     	long_name         vegetation temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   TWS                    	long_name         total water storage    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   U10                    	long_name         	10-m wind      units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   U10_DUST                   	long_name         10-m wind for dust model   units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   URBAN_AC                   	long_name         urban air conditioning flux    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   
URBAN_HEAT                     	long_name         urban heating flux     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   VCMX25T                    	long_name         canopy profile of vcmax25      units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   VENTILATION                    	long_name         ,sensible heat flux from building ventilation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   VOLR                   	long_name         !river channel total water storage      units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   VOLRMCH                    	long_name         (river channel main channel water storage   units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   VOLUMETRIC_STREAMFLOW                      	long_name         $volumetric streamflow from hillslope   units         m3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg             X   VPD2M                      	long_name         2m vapor pressure deficit      units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   VPD_CAN                    	long_name         canopy vapor pressure deficit      units         kPa    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   Vcmx25Z                    	long_name         1canopy profile of vcmax25 predicted by LUNA model      units         	umol/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   	WASTEHEAT                      	long_name         Csensible heat flux from heating/cooling sources of urban waste heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   WBT                    	long_name         2 m Stull Wet Bulb     units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   WBT_R                      	long_name         Rural 2 m Stull Wet Bulb   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   WBT_U                      	long_name         Urban 2 m Stull Wet Bulb   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   WIND                   	long_name         #atmospheric wind velocity magnitude    units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   ZBOT                   	long_name         atmospheric reference height   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   ZWT                    	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   	ZWT_PERCH                      	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   PCT_CFT                       	long_name         #% of each crop on the crop landunit    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   PCT_GLC_MEC                       	long_name         5% of each GLC elevation class on the glacier landunit      units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       (     �   SMP                       	long_name         Asoil matric potential (natural vegetated and crop landunits only)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d     �   TSOI                      	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d        TSOI_ICE                      	long_name         %soil temperature (ice landunits only)      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice       d     �   TLAKE            	             	long_name         lake temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       (     �   H2OSOI                        	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P        SOILICE                       	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     \   SOILLIQ                       	long_name         =soil liquid water (natural vegetated and crop landunits only)      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     �   PCT_LANDUNIT                      	long_name         % of each landunit on grid cell    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       $     �   PCT_NAT_PFT                       	long_name         =% of each PFT on the natural vegetation (i.e., soil) landunit      units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       <         VEGWP                         	long_name         Fvegetation water matric potential for sun/sha canopy,xyl,root segments     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   VEGWPLN                       	long_name         Kvegetation water matric potential for sun/sha canopy,xyl,root at local noon    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   VEGWPPD                       	long_name         Epredawn vegetation water matric potential for sun/sha canopy,xyl,root      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  C�A�B>@�=xČ?�           C� 1�E      �      Z�@|p     @~P     08/20/24        13:53:04        C�Ws+'��?��>),�-�[(            {@��A�I�A�I�>Q�5)$����)	�.��&0Z�>���>�|�    @6��A��A�k'Ak�>���B��B��C��pC��pC���>���B�w�    C��Bqn_CJ��Aw�BDaaC- 'AգOB��C��A�ARl�7�@�A��*`A�!A��A]T>=�>?#A���A1 -B
AV@_KL@�~A���@���G4|HG���G�i�G�H>S�Btu�A(>?��N���    @�t@�t{@��@Z{@Z{{@��BʵyBL��B���B�J>
D�=�E{@��{@��<mj�EJ�5?���*��(�i2*�vg*NL�*��(��+�>�+���+"pQ(�-�rB��gG��cB�';�E�;{᧴+�Z4�1�        7Q��    5�r�7hJ    83^	{@��1�5���8�����+{@��{@��{@��<���7��&            6��|8���8�8�{@��    8�    7m�6{@��8-�]{@��6cf6AEv5G�    7;��{@��43�5�/P7�'�7�'�B���B5b�    B��9A\�B#�4���3���>D�s6��l5�FA��	A��B#8:A��=�Y�=�F6|�:5OL�E#��7���=XF�Am�?_�M>��7���    8l�7��;C�,5BZ�Aa�Aa�{@�λ�#��#�C�֒{@��C�MvC���C�֒{@��>^�B��"EF�I@�͛C��C�v�C���>���C���C�AC��hC�39EP��?��?��i        B�<            ;|pC�wq>�mB��    ?��V?��V{@��?I+�A   ?/�r@.(�                                                à0yêhò��ì�â�~ÖЮÒ̘�x-�GQ��������H`��H��b,��@\��@\��@\��@\��@\��@\̾� ̾� ̾� ̾� ̾� C�AC��C��OC�e0C�
oC��C���C�ƚC��#C��WC�2C�/�C�I�C�d�C���C��+C��C��C�?~C�_�C���C���C�|UC�rxC�p�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��?e"?�0?t�?	�g?[�?�<>���>��[?�?��>�ɷ>�	�>��T>�X>��3>��3>��3>��3>��3>��3?)�9>˕�@'��@��A�F�B�Ao�k    =<n9?�	�                                        A&j�A�a�A��Bp�B#�B�B�5B�%CD�C![�C!�C*W�C�)C9��C]�C�O
C���C���C�SgCŪ1B�                                          @$                     Aᢌ        A���    @$     �Z�|�Y ��W���"�����ƽ�nƻ����Q��e�ƣ��e��;�