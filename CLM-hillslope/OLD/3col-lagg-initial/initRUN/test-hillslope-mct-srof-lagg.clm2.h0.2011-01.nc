CDF      
      lndgrid       gridcell      landunit      column        pft       levgrnd       levsoi        levurb        levmaxurbgrnd         levlak     
   numrad        levsno        ltype      	   nlevcan       nvegwcs       
nhillslope        max_columns_hillslope         	mxsowings         
mxharvests        natpft        cft       glc_nec    
   elevclas      string_length         scale_type_string_length       levdcmp       hist_interval         time          %   title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     Conventions       CF-1.0     history       created on 08/22/23 11:52:23   source        #Community Terrestrial Systems Model    hostname      cheyenne   username      marielj    version       unknown    revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       test-hillslope-mct-srof-lagg   Surface_dataset       Tsurfdata_1x1pt_US-MBP_hist_16pfts_Irrig_CMIP6_simyr2000_HAND_3_col_hillslope_lagg.nc   Initial_conditions_dataset        finidat_interp_dest.nc     #PFT_physiological_constants_dataset       clm50_params.c211112.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1      
   levgrnd                	long_name         coordinate ground levels   units         m         d     �   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m         P     L   levlak        	         	long_name         coordinate lake levels     units         m         (     �   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m              �   time               	long_name         time   units         days since 2011-01-01 00:00:00     calendar      noleap     bounds        time_bounds            ,   mcdate                 	long_name         current date (YYYYMMDD)            0   mcsec                  	long_name         current seconds of current date    units         s              4   mdcur                  	long_name         current day (from base day)            8   mscur                  	long_name         current seconds of current day             <   nstep                  	long_name         	time step              @   time_bounds                   	long_name         history time interval endpoints            D   date_written                            T   time_written                            d   lon                 	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��           �   lat                 	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��           �   area                	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��           �   landfrac                	long_name         land fraction      
_FillValue        {@��   missing_value         {@��           �   landmask                	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����           �   pftmask                 	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����           �   nbedrock                	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����           �   ZSOI                   	long_name         
soil depth     units         m      
_FillValue        {@��   missing_value         {@��   landunit_mask         nonurb        d     �   DZSOI                      	long_name         soil thickness     units         m      
_FillValue        {@��   missing_value         {@��   landunit_mask         nonurb        d     H   WATSAT                     	long_name         'saturated soil water content (porosity)    units         mm3/mm3    
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d     �   SUCSAT                     	long_name         saturated soil matric potential    units         mm     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d        BSW                    	long_name         #slope of soil water retention curve    units         unitless   
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d     t   HKSAT                      	long_name          saturated hydraulic conductivity   units         mm s-1     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d     �   ZLAKE         	             	long_name         lake layer node depth      units         m      
_FillValue        {@��   missing_value         {@��   landunit_mask         lake      (     <   DZLAKE        	             	long_name         lake layer thickness   units         m      
_FillValue        {@��   missing_value         {@��   landunit_mask         lake      (     d   PCT_SAND                   	long_name         percent sand   units         percent    
_FillValue        {@��   missing_value         {@��      P     �   PCT_CLAY                   	long_name         percent clay   units         percent    
_FillValue        {@��   missing_value         {@��      P     �   ATM_TOPO                   	long_name         atmospheric surface height     units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   BCDEP                      	long_name         -total BC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   BTRANMN                    	long_name         *daily minimum of transpiration beta factor     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            |   DSL                    	long_name         dry surface layer thickness    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   DSTDEP                     	long_name         /total dust deposition (dry+wet) from atmosphere    units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   DSTFLXT                    	long_name         total surface dust emission    units         kg/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	EFLXBUILD                      	long_name         Cbuilding heat flux from change in interior building air temperature    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   EFLX_DYNBAL                    	long_name         0dynamic land cover change conversion energy flux   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   EFLX_GRND_LAKE                     	long_name         Bnet heat flux into lake/snow surface, excluding light transmission     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   EFLX_LH_TOT                    	long_name         !total latent heat flux [+ to atm]      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   EFLX_LH_TOT_R                      	long_name         Rural total evaporation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   ELAI                   	long_name         !exposed one-sided leaf area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   ERRH2O                     	long_name         total water conservation error     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	ERRH2OSNO                      	long_name         &imbalance in snow depth (liquid water)     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   ERRSEB                     	long_name         !surface energy conservation error      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   ERRSOI                     	long_name         #soil/lake energy conservation error    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   ERRSOL                     	long_name         "solar radiation conservation error     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   ESAI                   	long_name         !exposed one-sided stem area index      units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FCEV                   	long_name         canopy evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FCOV                   	long_name         fractional impermeable area    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   FCTR                   	long_name         canopy transpiration   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FGEV                   	long_name         ground evaporation     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FGR                    	long_name         Oheat flux into soil/snow including snow melt and lake / snow light transmission    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FGR12                      	long_name         %heat flux between soil layers 1 and 2      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FH2OSFC                    	long_name         +fraction of ground covered by surface water    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FIRA                   	long_name         !net infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FIRA_R                     	long_name         'Rural net infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FIRE                   	long_name         %emitted infrared (longwave) radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FIRE_R                     	long_name         +Rural emitted infrared (longwave) radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FLDS                   	long_name         Matmospheric longwave radiation (downscaled for glacier and hillslope columns)      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FPSN                   	long_name         photosynthesis     units         umol m-2 s-1   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSA                    	long_name         absorbed solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSAT                   	long_name         +fractional area with water table at surface    units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   FSDS                   	long_name         Satmospheric incident solar radiation (downscaled for glacier and hillslope columns)    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSND                     	long_name         #direct nir incident solar radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   FSDSNDLN                   	long_name         1direct nir incident solar radiation at local noon      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                FSDSNI                     	long_name         $diffuse nir incident solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               FSDSVD                     	long_name         #direct vis incident solar radiation    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               FSDSVDLN                   	long_name         1direct vis incident solar radiation at local noon      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               FSDSVI                     	long_name         $diffuse vis incident solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               FSDSVILN                   	long_name         2diffuse vis incident solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               FSDS_from_atm                      	long_name         Oatmospheric incident solar radiation received from atmosphere (pre-downscaling)    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               FSH                    	long_name         Ssensible heat not including correction for land use change and rain/snow conversion    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               FSH_G                      	long_name         sensible heat from ground      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                FSH_PRECIP_CONVERSION                      	long_name         ;Sensible heat flux from conversion of rain/snow atm forcing    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   FSH_R                      	long_name         Rural sensible heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   FSH_RUNOFF_ICE_TO_LIQ                      	long_name         Dsensible heat flux generated from conversion of ice runoff to liquid   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   FSH_TO_COUPLER                     	long_name         �sensible heat sent to coupler (includes corrections for land use change, rain/snow conversion and conversion of ice runoff to liquid)      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   FSH_V                      	long_name         sensible heat from veg     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   FSM                    	long_name         snow melt heat flux    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   FSNO                   	long_name         "fraction of ground covered by snow     units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   FSNO_EFF                   	long_name         ,effective fraction of ground covered by snow   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   FSR                    	long_name         reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   FSRND                      	long_name         $direct nir reflected solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   FSRNDLN                    	long_name         2direct nir reflected solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   FSRNI                      	long_name         %diffuse nir reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   FSRVD                      	long_name         $direct vis reflected solar radiation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   FSRVDLN                    	long_name         2direct vis reflected solar radiation at local noon     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   FSRVI                      	long_name         %diffuse vis reflected solar radiation      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   GSSHA                      	long_name          shaded leaf stomatal conductance   units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   GSSHALN                    	long_name         .shaded leaf stomatal conductance at local noon     units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   GSSUN                      	long_name          sunlit leaf stomatal conductance   units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   GSSUNLN                    	long_name         .sunlit leaf stomatal conductance at local noon     units         umol H20/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   H2OCAN                     	long_name         intercepted water      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   H2OSFC                     	long_name         surface water depth    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   H2OSNO                     	long_name         snow depth (liquid water)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            x   
H2OSNO_TOP                     	long_name         mass of snow in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   H2OSOI                        	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     �   HEAT_CONTENT1                      	long_name         #initial gridcell total heat content    units         J/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   HEAT_FROM_AC                   	long_name         Lsensible heat flux put into canyon due to heat removed from air conditioning   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   HIA                    	long_name         2 m NWS Heat Index     units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   HIA_R                      	long_name         Rural 2 m NWS Heat Index   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   HIA_U                      	long_name         Urban 2 m NWS Heat Index   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   HUMIDEX                    	long_name         2 m Humidex    units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	HUMIDEX_R                      	long_name         Rural 2 m Humidex      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	HUMIDEX_U                      	long_name         Urban 2 m Humidex      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   ICE_CONTENT1                   	long_name         "initial gridcell total ice content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   IWUELN                     	long_name         )local noon intrinsic water use efficiency      units         umolCO2/molH2O     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   JMX25T                     	long_name         canopy profile of jmax     units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   Jmx25Z                     	long_name         Bmaximum rate of electron transport at 25 Celcius for canopy layers     units         umol electrons/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   LAISHA                     	long_name          shaded projected leaf area index   units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                LAISUN                     	long_name          sunlit projected leaf area index   units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               LAKEICEFRAC_SURF                   	long_name         $surface lake layer ice mass fraction   units         unitless   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               LAKEICETHICK                   	long_name         @thickness of lake ice (including physical expansion on freezing)   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               LIQCAN                     	long_name         intercepted liquid water   units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               LIQUID_CONTENT1                    	long_name         "initial gridcell total liq content     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               LNC                    	long_name         leaf N concentration   units         gN leaf/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               MEG_acetaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               MEG_acetic_acid                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                MEG_acetone                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   MEG_carene_3                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   MEG_ethanol                    	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   MEG_formaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   MEG_isoprene                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   MEG_methanol                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   MEG_pinene_a                   	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   MEG_thujene_a                      	long_name         
MEGAN flux     units         	kg/m2/sec      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   OCDEP                      	long_name         -total OC deposition (dry+wet) from atmosphere      units         kg/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   PARVEGLN                   	long_name         (absorbed par by vegetation at local noon   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   PBOT                   	long_name         Natmospheric pressure at surface (downscaled for glacier and hillslope columns)     units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   PCO2                   	long_name         #atmospheric partial pressure of CO2    units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   PCT_CFT                       	long_name         #% of each crop on the crop landunit    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   PCT_GLC_MEC                       	long_name         5% of each GLC elevation class on the glacier landunit      units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       (     \   PCT_LANDUNIT                      	long_name         % of each landunit on grid cell    units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       $     �   PCT_NAT_PFT                       	long_name         =% of each PFT on the natural vegetation (i.e., soil) landunit      units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       <     �   Q2M                    	long_name         2m specific humidity   units         kg/kg      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QBOT                   	long_name         Hatmospheric specific humidity (downscaled to columns in glacier regions)   units         kg/kg      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QDRAI                      	long_name         sub-surface drainage   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QDRAI_PERCH                    	long_name         perched wt drainage    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QDRAI_XS                   	long_name         saturation excess drainage     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLOOD                     	long_name         runoff from river flooding     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_EVAP_TOT                      	long_name         -qflx_evap_soi + qflx_evap_can + qflx_tran_veg      units         
kg m-2 s-1     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QFLX_ICE_DYNBAL                    	long_name         4ice dynamic land cover change conversion runoff flux   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                QFLX_LIQDEW_TO_TOP_LAYER                   	long_name         >rate of liquid water deposited on top soil or snow layer (dew)     units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QFLX_LIQEVAP_FROM_TOP_LAYER                    	long_name         ;rate of liquid water evaporated from top soil or snow layer    units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QFLX_LIQ_DYNBAL                    	long_name         4liq dynamic land cover change conversion runoff flux   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QFLX_SNOW_DRAIN                    	long_name         drainage from snow pack    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QFLX_SNOW_DRAIN_ICE                    	long_name         1drainage from snow pack melt (ice landunits only)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice               QFLX_SOLIDDEW_TO_TOP_LAYER                     	long_name         ?rate of solid water deposited on top soil or snow layer (frost)    units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QFLX_SOLIDEVAP_FROM_TOP_LAYER                      	long_name         zrate of ice evaporated from top soil or snow layer (sublimation) (also includes bare ice sublimation from glacier columns)     units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               QH2OSFC                    	long_name         surface water runoff   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                QHR                    	long_name         hydraulic redistribution   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            $   QICE                   	long_name         ice growth/melt    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            (   QICE_FRZ                   	long_name         
ice growth     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            ,   	QICE_MELT                      	long_name         ice melt   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            0   QINFL                      	long_name         infiltration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   QINTR                      	long_name         interception   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   QIRRIG_FROM_GW_CONFINED                    	long_name         3water added through confined groundwater irrigation    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   QIRRIG_FROM_GW_UNCONFINED                      	long_name         5water added through unconfined groundwater irrigation      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   QIRRIG_FROM_SURFACE                    	long_name         ,water added through surface water irrigation   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   QLATFLOWOUT                    	long_name         hillcol lateral outflow    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg             H   QOVER                      	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   QRGWL                      	long_name         isurface runoff at glaciers (liquid only), wetlands, lakes; also includes melted ice runoff from QSNWCPICE      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   QRUNOFF                    	long_name         @total liquid runoff not including correction for land use change   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   QRUNOFF_ICE                    	long_name         Btotal liquid runoff not incl corret for LULCC (ice landunits only)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            X   QRUNOFF_ICE_TO_COUPLER                     	long_name         Ktotal ice runoff sent to coupler (includes corrections for land use change)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   QRUNOFF_TO_COUPLER                     	long_name         Ntotal liquid runoff sent to coupler (includes corrections for land use change)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   	QSNOCPLIQ                      	long_name         Rexcess liquid h2o due to snow capping not including correction for land use change     units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   QSNOEVAP                   	long_name         Nevaporation from snow (only when snl<0, otherwise it is equal to qflx_ev_soil)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   QSNOFRZ                    	long_name         $column-integrated snow freezing rate   units         kg/m2/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            l   QSNOFRZ_ICE                    	long_name         9column-integrated snow freezing rate (ice landunits only)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            p   QSNOMELT                   	long_name         snow melt rate     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   QSNOMELT_ICE                   	long_name         snow melt (ice landunits only)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            x   
QSNOUNLOAD                     	long_name         canopy snow unloading      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            |   QSNO_TEMPUNLOAD                    	long_name         canopy snow temp unloading     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QSNO_WINDUNLOAD                    	long_name         canopy snow wind unloading     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	QSNWCPICE                      	long_name         Qexcess solid h2o due to snow capping not including correction for land use change      units         mm H2O/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QSOIL                      	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   	QSOIL_ICE                      	long_name         'Ground evaporation (ice landunits only)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice            �   QVEGE                      	long_name         canopy evaporation     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   QVEGT                      	long_name         canopy transpiration   units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   RAIN                   	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   RAIN_FROM_ATM                      	long_name         >atmospheric rain received from atmosphere (pre-repartitioning)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   RH2M                   	long_name         2m relative humidity   units         %      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   RSSHA                      	long_name         shaded leaf stomatal resistance    units         s/m    cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   RSSUN                      	long_name         sunlit leaf stomatal resistance    units         s/m    cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   SABG                   	long_name         solar rad absorbed by ground   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SABG_PEN                   	long_name         2Rural solar rad penetrating top soil or snow layer     units         watt/m^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SABV                   	long_name         solar rad absorbed by veg      units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   SMP                       	long_name         Asoil matric potential (natural vegetated and crop landunits only)      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d     �   SNOBCMCL                   	long_name         mass of BC in snow column      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                SNOBCMSL                   	long_name         mass of BC in top snow layer   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   SNOCAN                     	long_name         intercepted snow   units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   	SNODSTMCL                      	long_name         mass of dust in snow column    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   	SNODSTMSL                      	long_name         mass of dust in top snow layer     units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   SNOFSRND                   	long_name         .direct nir reflected solar radiation from snow     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   SNOFSRNI                   	long_name         /diffuse nir reflected solar radiation from snow    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   SNOFSRVD                   	long_name         .direct vis reflected solar radiation from snow     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   SNOFSRVI                   	long_name         /diffuse vis reflected solar radiation from snow    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   	SNOINTABS                      	long_name         8Fraction of incoming solar absorbed by lower snow layers   units         -      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   SNOOCMCL                   	long_name         mass of OC in snow column      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   SNOOCMSL                   	long_name         mass of OC in top snow layer   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   	SNOTXMASS                      	long_name         ksnow temperature times layer mass, layer sum; to get mass-weighted temperature, divide by (SNOWICE+SNOWLIQ)    units         K kg/m2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            P   SNOW                   	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            T   SNOWDP                     	long_name         gridcell mean snow height      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            X   SNOWICE                    	long_name         snow ice   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            \   SNOWLIQ                    	long_name         snow liquid water      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            `   
SNOW_DEPTH                     	long_name          snow height of snow covered area   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            d   SNOW_FROM_ATM                      	long_name         >atmospheric snow received from atmosphere (pre-repartitioning)     units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            h   SNOW_PERSISTENCE                   	long_name         BLength of time of continuous snow cover (nat. veg. landunits only)     units         seconds    cell_methods      time: instantaneous    
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg             l   
SNOW_SINKS                     	long_name         snow sinks (liquid water)      units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            p   SNOW_SOURCES                   	long_name         snow sources (liquid water)    units         mm/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            t   SOILICE                       	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     x   SOILLIQ                       	long_name         =soil liquid water (natural vegetated and crop landunits only)      units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       P     �   	SOILRESIS                      	long_name         soil resistance to evaporation     units         s/m    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown               SOILWATER_10CM                     	long_name         @soil liquid water + ice in top 10cm of soil (veg landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg               SWBGT                      	long_name         !2 m Simplified Wetbulb Globe Temp      units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                SWBGT_R                    	long_name         'Rural 2 m Simplified Wetbulb Globe Temp    units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            $   SWBGT_U                    	long_name         'Urban 2 m Simplified Wetbulb Globe Temp    units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            (   TAUX                   	long_name         zonal surface stress   units         kg/m/s^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            ,   TAUY                   	long_name         meridional surface stress      units         kg/m/s^2   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            0   TBOT                   	long_name         Jatmospheric air temperature (downscaled for glacier and hillslope columns)     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            4   TBUILD                     	long_name         'internal urban building air temperature    units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            8   TG                     	long_name         ground temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            <   TH2OSFC                    	long_name         surface water temperature      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            @   THBOT                      	long_name         Tatmospheric air potential temperature (downscaled for glacier and hillslope columns)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            D   TKE1                   	long_name         (top lake level eddy thermal conductivity   units         W/(mK)     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            H   TLAI                   	long_name         total projected leaf area index    units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            L   TLAKE            	             	long_name         lake temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       (     P   
TOTSOILICE                     	long_name         /vertically summed soil cie (veg landunits only)    units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            x   
TOTSOILLIQ                     	long_name         8vertically summed soil liquid water (veg landunits only)   units         kg/m2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            |   TPU25T                     	long_name         canopy profile of tpu      units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TREFMNAV                   	long_name         (daily minimum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TREFMXAV                   	long_name         (daily maximum of average 2-m temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TSA                    	long_name         2m air temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TSAI                   	long_name         total projected stem area index    units         m^2/m^2    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TSKIN                      	long_name         skin temperature   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            �   TSL                    	long_name         Rtemperature of near-surface soil layer (natural vegetated and crop landunits only)     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg            �   TSOI                      	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg       d     �   	TSOI_10CM                      	long_name         $soil temperature in top 10cm of soil   units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown                 TSOI_ICE                      	long_name         %soil temperature (ice landunits only)      units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         ice       d         TV                     	long_name         vegetation temperature     units         K      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             h   TWS                    	long_name         total water storage    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             l   U10                    	long_name         	10-m wind      units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             p   U10_DUST                   	long_name         10-m wind for dust model   units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             t   URBAN_AC                   	long_name         urban air conditioning flux    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             x   
URBAN_HEAT                     	long_name         urban heating flux     units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             |   VCMX25T                    	long_name         canopy profile of vcmax25      units         	umol/m2/s      cell_methods      time: minimum      
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VEGWP                         	long_name         Fvegetation water matric potential for sun/sha canopy,xyl,root segments     units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VEGWPLN                       	long_name         Kvegetation water matric potential for sun/sha canopy,xyl,root at local noon    units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VEGWPPD                       	long_name         Epredawn vegetation water matric potential for sun/sha canopy,xyl,root      units         mm     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VENTILATION                    	long_name         ,sensible heat flux from building ventilation   units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VOLR                   	long_name         !river channel total water storage      units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VOLRMCH                    	long_name         (river channel main channel water storage   units         m3     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VOLUMETRIC_STREAMFLOW                      	long_name         $volumetric streamflow from hillslope   units         m3/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         natveg              �   VPD2M                      	long_name         2m vapor pressure deficit      units         Pa     cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   VPD_CAN                    	long_name         canopy vapor pressure deficit      units         kPa    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   Vcmx25Z                    	long_name         1canopy profile of vcmax25 predicted by LUNA model      units         	umol/m2/s      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   	WASTEHEAT                      	long_name         Csensible heat flux from heating/cooling sources of urban waste heat    units         W/m^2      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   WBT                    	long_name         2 m Stull Wet Bulb     units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   WBT_R                      	long_name         Rural 2 m Stull Wet Bulb   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   WBT_U                      	long_name         Urban 2 m Stull Wet Bulb   units         C      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   WIND                   	long_name         #atmospheric wind velocity magnitude    units         m/s    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   ZBOT                   	long_name         atmospheric reference height   units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             �   ZWT                    	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg             �   	ZWT_PERCH                      	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         veg             �<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�=L��?��@ff@�33A��AI��A���A���B	L�B3�?�  C�A�B>@�B��?�           <#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<��
=#�
=u=��
=�>#�
>L��>u>�\)>��
>�Q�>���>�G�?
=q?#�
?=p�?W
=?p��?��?��@��@���@�S{A2=qAq��?m�h?l1?ix�?e�T?`Ĝ?Y��?Tz�?Tz�?Tz�?Tz�?(�k?(�k>��>��>��>��>��>��>��>��>��>��>��>��>��A!��A!��A!��A!��A!��A!��A!�A�AA�C��dC��D2��D2��D2��D2��D2��D2��D2��D2��D2��D2��D2��D2��D2��@8�9@\j@���@���@�&�A"=qA@  A@  A@  A@  @��@��@N��@N��@N��@N��@N��@N��@N��@N��@N��@N��@N��@N��@N��>�~k>��2>k �>C �>	��=e�$:���:���:���:���:���:���:���:���:���:���:���:���:���:���:���:���:���:���:���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��@   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   @   A�  2�y             п�UUUUUU@?      08/22/23        11:52:23        CŢ3+�?m��@ס�+tM�            {@��@@    '"�'�T&�&#
r/2�#JcO?+Q>�c�=�q�%N�Y?�Q�_;^�`&�6A�m@��@��Cz��Cz��Cu��    A'��=�q�A��o@�oYA�h�@���@;�Apڡ@�A�4�A��oA���A\�    A���    A���@tL�    ?Y|?Y|@���?�h�@�@�?��~?�7JA��@ZRFtUSF�7�F��?F�]d?�\x6�[%A�4?Ɏ�?�9?�>>�~x?��?%X?#tD?X�?&��?g�?�y? �f?��>���>���>��#>���>�e<>��>��>��N��P    ��t���t�{@����6���6�{@��C��{@��B
XBp          {@��{@��    Ec��?���                                        ,��OA@ �G��[B�X                                                B�                                                              BMh                @��        :~c�:y�����U            5E�3                    {@��0�i�5�    ����{@��{@��{@�δ�16w��            0��        ���U{@��    ���U    5��    {@��    {@��6E�     6E�     5�y{@��47�����        B�a        @h�@!2A@ڡ��\�K�E�>�����8��F��7���s\�9�f�����_��0��$��N���<��"���hġ]��2���2���2��̾� ̾� ̾� ̾� ̾� 5��n4 �&?�\x6;�4a�@3s�@L8@(�@���>J�7q�u5�ݍE��7\ԝ>�A�4    >x�7\ԝJx��4�_�7Pn�A�aA�A�A�(BE?BN�BU�!BZ�@�1�                                                <vݘ?<��@C�AE�A���B/^B���C$�C(�C"]�C5
WCL-C'C(ESCP�kC��C�!C�FC��8D
��C�M�BC�
��#���#�{@�λ�!b��!bC�8�{@��C��#C���C�8�{@��;2��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Ck�
Ec��@@�KCz�/C�W�C�[�?,$2C���C�nOC�nOC���C���C�{zC�8�C�j�C�aWC���C�ΥC���C�7HC�~=C���C�ĀC�գC���C��nC�ɱC��C���C��JC��ZC���C��C��iC���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��!Et?��?U�        A�AƔ6�Ɣ6�Ɣ6���/�ƞ��ƞ��ƞ���ζ����������ě2?            ���A�0w=j��A�      ������{@��?
g�A   @�Z�    