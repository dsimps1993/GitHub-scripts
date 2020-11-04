##########################################################
### iDat processing and age pred      04/11/20 DJS
##########################################################


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bigmelon",'wateRmelon', 'gdsfmt')

library(bigmelon)
gfile <- iadd2('GSE61461_RAW/', gds = 'GSE61461.gds', chunksize = 100)

index.gdsn(gfile, 'betas')
class(index.gdsn(gfile, 'betas'))
index.gdsn(gfile, 'fData/TargetID')

gfile <- openfn.gds("GSE61461.gds")

###extracting the horvy cpg stuff
##these are the age cpgs (Nelly I can send them to you later, but you can also download them from Horv's 2013 paper)
load("C:/Users/s1684303/Google Drive/University of Edinburgh PhD/18 CpG Project/R files/coef.rda")


#agep function applies the Horvath age predictor to your data.
agep(betas(gfile), coef)

# THis metadata object was one I downloaded (see code below) for this dataset. You might have one come with your data or you might have to make your own. Essentially all it needs to have is the samples as rownames (matching order of the samples in the array), and the columns will be data like chronological age, gender etc, and the age predicitons you'll be adding later.

library(GEOquery)
GSE61461.meta <- getGEO("GSE61461", GSEMatrix = TRUE)

GSE61461.metaDat <- GSE61461.meta[["GSE61461_series_matrix.txt.gz"]]@phenoData@data
write.csv(GSE61461.metaDat, "GSE61461.metaDat.csv")

### reading in new dat with extra columns for plotting later
GSE61461.metaDat2 <- read.csv("GSE61461.metaDat2.csv", header = T, row.names = 1)




##adding in eAge to meta 

GSE61461.metaDat2$Horvath <- as.numeric(agep(betas(gfile), coef))

####Horvath Skin and Blood clock #############

###it seems Im struggling to use the gfile for this age pred, because the clock requires the rownames to be the first column of the dataframe
###so, I'm gonna save as a csv and load it again.

##Square brackets with comma here converts it to a dataframe ish object
write.csv(betas(gfile)[,], "GSE614161_Fullarray.csv")

#This closes the gfile, its hard to explain exactly why theres an open and closed gfile thing, but its all good because you have the matrix saved as .csv and can load them up.
closefn.gds(gfile)

dat0 <- read.csv("GSE614161_Fullarray.csv", header = T)


#R functions for transforming age
adult.age1=20
trafo= function(x,adult.age=adult.age1)
{x=(x+1)/(1+adult.age); y=ifelse(x<=1,
                                 log( x),x-1);y }
anti.trafo=
  function(x,adult.age=adult.age1) {
    ifelse(x<0, (1+adult.age)*exp(x)-1,
           (1+adult.age)*x+adult.age) }
datClock=read.csv("HorvSkin&Bld_clock.csv")


###na turned to zero version, if you want I can send you code where you use mean CpG methylation to impute instead
#dat0[is.na(dat0)] <- 0


selectCpGsClock=is.element(as.character(dat0[,1]),
                           as.character(datClock[-1,1]))
datMethClock0=data.frame(t(dat0[selectCpGsClock ,-1] ))

colnames(datMethClock0)=
  as.character(dat0[selectCpGsClock ,1] )
# Reality check: the following output should only contain numeric values.
# Further, the column names should be CpG identifiers (cg numbers).
datMethClock0[1:5,1:5]
datMethClock= data.frame(datMethClock0[
  as.character(datClock[-1,1])])
# The number of rows should equal the number of samples (Illumina arrays)
dim(datMethClock)
#Output DNAm age estimator for the skin & blood clock
DNAmAgeSkinClock=as.numeric(anti.trafo(
  datClock$Coef[1]+as.matrix(datMethClock
  )%*% as.numeric(datClock$Coef[-1])))


##adding in age to meta
GSE61461.metaDat2$SkinClock <- DNAmAgeSkinClock





#### PhenoAge clock #####################

##not memory conscious here, but reading in the object again with the CpGs as rownames.
datMeth <- t(read.csv("GSE614161_Fullarray.csv", header = T, row.names = 1))
datMeth <- as.data.frame(datMeth, stringsAsFactors = FALSE)



###reading in phenoageCpgs and extracting em
DNAmPhenoAge_CpGs <- read.table(file.choose(), header = T)


datMeth2 <- datMeth[,colnames(datMeth) %in% DNAmPhenoAge_CpGs$ID_REF]


###this is just checking which are NAs
for(i in colnames(datMeth2)){
  print(i)
  p = which(is.na(datMeth2[,i]) == TRUE)
  print(p)
  for(n in p){
    datMeth2[n,i] <- mean(as.numeric(datMeth2[,i]), na.rm = T)
    print(datMeth2[n,i])
  }
}

datMeth <- as.data.frame(datMeth2)




DNAmPhenoAge=60.66387252+
  (0.122818203*as.numeric(datMeth$cg00079056))+
  (0.030577331*as.numeric(datMeth$cg00083937))+
  (4.02670329*as.numeric(datMeth$cg00113951))+
  (-1.375221101*as.numeric(datMeth$cg00168942))+
  (-11.88629118*as.numeric(datMeth$cg00194146))+
  (-12.15471329*as.numeric(datMeth$cg00230271))+
  (13.39937728*as.numeric(datMeth$cg00261781))+
  (1.908695219*as.numeric(datMeth$cg00297600))+
  (1.855444605*as.numeric(datMeth$cg00335286))+
  (0.534125992*as.numeric(datMeth$cg00338702))+
  (2.072016695*as.numeric(datMeth$cg00350702))+
  (-4.303054603*as.numeric(datMeth$cg00410898))+
  (1.844120976*as.numeric(datMeth$cg00412772))+
  (4.016808807*as.numeric(datMeth$cg00412805))+
  (13.58737822*as.numeric(datMeth$cg00462994))+
  (0.002679625*as.numeric(datMeth$cg00503840))+
  (-4.103268356*as.numeric(datMeth$cg00515905))+
  (-1.678163605*as.numeric(datMeth$cg00582628))+
  (3.757669266*as.numeric(datMeth$cg00687674))+
  (-1.411402743*as.numeric(datMeth$cg00744433))+
  (-8.58936172*as.numeric(datMeth$cg00845900))+
  (-0.233766711*as.numeric(datMeth$cg00862290))+
  (33.47697933*as.numeric(datMeth$cg00943950))+
  (3.41792265*as.numeric(datMeth$cg00955230))+
  (-0.326099256*as.numeric(datMeth$cg01056568))+
  (-0.226611181*as.numeric(datMeth$cg01114088))+
  (-0.499180939*as.numeric(datMeth$cg01128603))+
  (-22.91466682*as.numeric(datMeth$cg01131735))+
  (14.37377079*as.numeric(datMeth$cg01137065))+
  (-10.27522152*as.numeric(datMeth$cg01211097))+
  (-8.838557811*as.numeric(datMeth$cg01221637))+
  (2.348957061*as.numeric(datMeth$cg01252496))+
  (19.18438139*as.numeric(datMeth$cg01254459))+
  (1.047474154*as.numeric(datMeth$cg01261503))+
  (1.211111994*as.numeric(datMeth$cg01335367))+
  (4.009589827*as.numeric(datMeth$cg01400401))+
  (-8.741688217*as.numeric(datMeth$cg01441777))+
  (-11.65425604*as.numeric(datMeth$cg01450842))+
  (-0.294798529*as.numeric(datMeth$cg01459453))+
  (-0.955659544*as.numeric(datMeth$cg01511567))+
  (7.247589771*as.numeric(datMeth$cg01519742))+
  (-3.806275602*as.numeric(datMeth$cg01623187))+
  (-5.640661125*as.numeric(datMeth$cg01626227))+
  (12.43450304*as.numeric(datMeth$cg01651821))+
  (0.481129756*as.numeric(datMeth$cg01918706))+
  (-21.78792705*as.numeric(datMeth$cg01930621))+
  (-1.552392375*as.numeric(datMeth$cg01946401))+
  (7.246984432*as.numeric(datMeth$cg02016419))+
  (-3.886409522*as.numeric(datMeth$cg02071305))+
  (-6.962548699*as.numeric(datMeth$cg02151301))+
  (6.199036211*as.numeric(datMeth$cg02154074))+
  (1.326475211*as.numeric(datMeth$cg02197293))+
  (-7.785341156*as.numeric(datMeth$cg02228185))+
  (10.5122237*as.numeric(datMeth$cg02229946))+
  (-5.609942057*as.numeric(datMeth$cg02309431))+
  (1.531112896*as.numeric(datMeth$cg02480835))+
  (25.58514041*as.numeric(datMeth$cg02503970))+
  (-7.435453213*as.numeric(datMeth$cg02631957))+
  (4.59476785*as.numeric(datMeth$cg02735486))+
  (-24.11715679*as.numeric(datMeth$cg02802055))+
  (8.057451538*as.numeric(datMeth$cg02976574))+
  (-2.202237705*as.numeric(datMeth$cg03007010))+
  (9.309963576*as.numeric(datMeth$cg03112869))+
  (-11.96946328*as.numeric(datMeth$cg03172991))+
  (-3.301812331*as.numeric(datMeth$cg03258472))+
  (13.57369508*as.numeric(datMeth$cg03340261))+
  (4.630001389*as.numeric(datMeth$cg03387497))+
  (-15.95141109*as.numeric(datMeth$cg03535648))+
  (7.763007186*as.numeric(datMeth$cg03565081))+
  (1.943750756*as.numeric(datMeth$cg03623878))+
  (2.057604969*as.numeric(datMeth$cg03703325))+
  (6.080591714*as.numeric(datMeth$cg03724882))+
  (4.145197074*as.numeric(datMeth$cg03819692))+
  (-19.41008226*as.numeric(datMeth$cg03929796))+
  (19.45230678*as.numeric(datMeth$cg03977782))+
  (-7.60434303*as.numeric(datMeth$cg03991512))+
  (-8.441092863*as.numeric(datMeth$cg04007936))+
  (-7.924579894*as.numeric(datMeth$cg04014889))+
  (31.62165709*as.numeric(datMeth$cg04084157))+
  (2.383062953*as.numeric(datMeth$cg04087608))+
  (4.926923884*as.numeric(datMeth$cg04169469))+
  (0.164994457*as.numeric(datMeth$cg04333463))+
  (0.255790471*as.numeric(datMeth$cg04359302))+
  (-5.338257391*as.numeric(datMeth$cg04416752))+
  (-14.44777739*as.numeric(datMeth$cg04424621))+
  (-4.215631426*as.numeric(datMeth$cg04480914))+
  (11.07259285*as.numeric(datMeth$cg04528819))+
  (-4.847499211*as.numeric(datMeth$cg04601137))+
  (0.068531735*as.numeric(datMeth$cg04616566))+
  (-5.452605586*as.numeric(datMeth$cg04718414))+
  (-16.11768623*as.numeric(datMeth$cg04736140))+
  (3.817194129*as.numeric(datMeth$cg04755031))+
  (0.530063706*as.numeric(datMeth$cg04818845))+
  (10.65458827*as.numeric(datMeth$cg04836038))+
  (-6.991225132*as.numeric(datMeth$cg05087948))+
  (-3.91160807*as.numeric(datMeth$cg05089968))+
  (-4.602562372*as.numeric(datMeth$cg05125838))+
  (-4.49429672*as.numeric(datMeth$cg05228408))+
  (9.458683059*as.numeric(datMeth$cg05270634))+
  (10.43369407*as.numeric(datMeth$cg05294243))+
  (-14.22870517*as.numeric(datMeth$cg05316065))+
  (0.965636264*as.numeric(datMeth$cg05422352))+
  (1.547194374*as.numeric(datMeth$cg05440289))+
  (3.103760739*as.numeric(datMeth$cg05441133))+
  (-3.966098185*as.numeric(datMeth$cg05442902))+
  (-5.15604937*as.numeric(datMeth$cg05473871))+
  (-16.73113119*as.numeric(datMeth$cg05492270))+
  (-0.164581823*as.numeric(datMeth$cg05501584))+
  (6.729463914*as.numeric(datMeth$cg05532892))+
  (-0.297182577*as.numeric(datMeth$cg05697249))+
  (-9.047212407*as.numeric(datMeth$cg05759269))+
  (-8.056495421*as.numeric(datMeth$cg05851163))+
  (1.202033611*as.numeric(datMeth$cg05898102))+
  (5.64268843*as.numeric(datMeth$cg06134964))+
  (3.413484589*as.numeric(datMeth$cg06144905))+
  (0.573872643*as.numeric(datMeth$cg06171242))+
  (-9.928634297*as.numeric(datMeth$cg06189653))+
  (1.251326977*as.numeric(datMeth$cg06295856))+
  (5.716566436*as.numeric(datMeth$cg06327515))+
  (-0.047794699*as.numeric(datMeth$cg06363129))+
  (23.37207957*as.numeric(datMeth$cg06493994))+
  (16.05623965*as.numeric(datMeth$cg06533629))+
  (1.89342443*as.numeric(datMeth$cg06637774))+
  (2.479689072*as.numeric(datMeth$cg06638451))+
  (-17.31209095*as.numeric(datMeth$cg06690548))+
  (10.86511161*as.numeric(datMeth$cg06908778))+
  (2.788051118*as.numeric(datMeth$cg06958034))+
  (-0.508782249*as.numeric(datMeth$cg06975499))+
  (2.264576345*as.numeric(datMeth$cg06994793))+
  (1.7787351*as.numeric(datMeth$cg07038400))+
  (-2.216579511*as.numeric(datMeth$cg07073964))+
  (-4.105965549*as.numeric(datMeth$cg07180649))+
  (-2.423609122*as.numeric(datMeth$cg07211259))+
  (4.684481441*as.numeric(datMeth$cg07236943))+
  (28.3980076*as.numeric(datMeth$cg07265300))+
  (6.581498688*as.numeric(datMeth$cg07484827))+
  (-5.794138096*as.numeric(datMeth$cg07494518))+
  (-7.183441358*as.numeric(datMeth$cg07654934))+
  (12.90897644*as.numeric(datMeth$cg07817698))+
  (21.31021579*as.numeric(datMeth$cg07850604))+
  (-9.267771972*as.numeric(datMeth$cg07929310))+
  (-6.589880588*as.numeric(datMeth$cg08035942))+
  (-11.51984975*as.numeric(datMeth$cg08067365))+
  (1.902730524*as.numeric(datMeth$cg08074477))+
  (-6.322778588*as.numeric(datMeth$cg08169325))+
  (0.012944061*as.numeric(datMeth$cg08212685))+
  (0.755043977*as.numeric(datMeth$cg08251399))+
  (-20.70742942*as.numeric(datMeth$cg08331960))+
  (-1.156840533*as.numeric(datMeth$cg08424423))+
  (3.702746231*as.numeric(datMeth$cg08475827))+
  (11.67344393*as.numeric(datMeth$cg08487374))+
  (-10.22704896*as.numeric(datMeth$cg08529529))+
  (2.696221083*as.numeric(datMeth$cg08586737))+
  (-3.06398296*as.numeric(datMeth$cg08587542))+
  (1.973711223*as.numeric(datMeth$cg08654655))+
  (4.112113212*as.numeric(datMeth$cg08668790))+
  (2.481527201*as.numeric(datMeth$cg08694544))+
  (1.548231761*as.numeric(datMeth$cg08872493))+
  (-16.3038063*as.numeric(datMeth$cg08896629))+
  (11.08380298*as.numeric(datMeth$cg08899632))+
  (2.536833226*as.numeric(datMeth$cg08900043))+
  (0.11706097*as.numeric(datMeth$cg09045681))+
  (0.684630536*as.numeric(datMeth$cg09096950))+
  (0.260322504*as.numeric(datMeth$cg09196959))+
  (-1.115123986*as.numeric(datMeth$cg09254939))+
  (-0.157799598*as.numeric(datMeth$cg09294589))+
  (3.548407354*as.numeric(datMeth$cg09304040))+
  (3.784159584*as.numeric(datMeth$cg09322949))+
  (-0.481760063*as.numeric(datMeth$cg09404633))+
  (-3.403045935*as.numeric(datMeth$cg09413557))+
  (0.04407787*as.numeric(datMeth$cg09434995))+
  (-7.741831686*as.numeric(datMeth$cg09480837))+
  (-19.0323049*as.numeric(datMeth$cg09548179))+
  (-18.61940178*as.numeric(datMeth$cg09556292))+
  (0.165360604*as.numeric(datMeth$cg09630437))+
  (-0.750651448*as.numeric(datMeth$cg09799873))+
  (-12.09921682*as.numeric(datMeth$cg09809672))+
  (2.879940739*as.numeric(datMeth$cg09851465))+
  (-33.4823461*as.numeric(datMeth$cg09892203))+
  (-5.67467637*as.numeric(datMeth$cg10052840))+
  (-10.1880066*as.numeric(datMeth$cg10158181))+
  (-0.41970649*as.numeric(datMeth$cg10202457))+
  (11.81083061*as.numeric(datMeth$cg10225525))+
  (5.941622921*as.numeric(datMeth$cg10523019))+
  (2.358758156*as.numeric(datMeth$cg10570177))+
  (-1.280172495*as.numeric(datMeth$cg10591174))+
  (-12.95489*as.numeric(datMeth$cg10636246))+
  (13.83949625*as.numeric(datMeth$cg10654016))+
  (-0.494137066*as.numeric(datMeth$cg10667970))+
  (-1.649159645*as.numeric(datMeth$cg10669058))+
  (-0.630165*as.numeric(datMeth$cg10795646))+
  (5.340752668*as.numeric(datMeth$cg10878896))+
  (-0.993401947*as.numeric(datMeth$cg10900550))+
  (-2.815298757*as.numeric(datMeth$cg10917602))+
  (10.79008118*as.numeric(datMeth$cg10922280))+
  (-10.94667271*as.numeric(datMeth$cg11177450))+
  (-8.702162118*as.numeric(datMeth$cg11233384))+
  (-16.243299*as.numeric(datMeth$cg11237115))+
  (-24.78579343*as.numeric(datMeth$cg11426590))+
  (-4.21297823*as.numeric(datMeth$cg11459714))+
  (-5.062657046*as.numeric(datMeth$cg11487705))+
  (1.961425559*as.numeric(datMeth$cg11490446))+
  (-15.57182248*as.numeric(datMeth$cg11600161))+
  (0.45590634*as.numeric(datMeth$cg11618577))+
  (6.461478745*as.numeric(datMeth$cg11631518))+
  (-5.210109944*as.numeric(datMeth$cg11833861))+
  (-1.375881082*as.numeric(datMeth$cg11896923))+
  (5.525340985*as.numeric(datMeth$cg11903057))+
  (2.227238015*as.numeric(datMeth$cg12145907))+
  (2.336851094*as.numeric(datMeth$cg12177001))+
  (-1.546553159*as.numeric(datMeth$cg12188560))+
  (5.749756441*as.numeric(datMeth$cg12238343))+
  (-20.75342001*as.numeric(datMeth$cg12247247))+
  (-8.80006302*as.numeric(datMeth$cg12261786))+
  (12.6784476*as.numeric(datMeth$cg12265604))+
  (10.4580988*as.numeric(datMeth$cg12269343))+
  (0.205379773*as.numeric(datMeth$cg12289045))+
  (0.299599643*as.numeric(datMeth$cg12324144))+
  (0.01052681*as.numeric(datMeth$cg12373771))+
  (1.509543012*as.numeric(datMeth$cg12402251))+
  (-3.994735863*as.numeric(datMeth$cg12473775))+
  (36.7881845*as.numeric(datMeth$cg12743894))+
  (1.182049823*as.numeric(datMeth$cg12813792))+
  (0.655304718*as.numeric(datMeth$cg12864235))+
  (-35.90008808*as.numeric(datMeth$cg12985418))+
  (-0.399519083*as.numeric(datMeth$cg12991365))+
  (-8.160794723*as.numeric(datMeth$cg13042288))+
  (-0.497800027*as.numeric(datMeth$cg13119609))+
  (-8.761639169*as.numeric(datMeth$cg13120519))+
  (-13.68132768*as.numeric(datMeth$cg13218906))+
  (1.935417196*as.numeric(datMeth$cg13258700))+
  (4.457073326*as.numeric(datMeth$cg13296371))+
  (2.583560626*as.numeric(datMeth$cg13307384))+
  (3.552242836*as.numeric(datMeth$cg13323474))+
  (-3.960753228*as.numeric(datMeth$cg13351161))+
  (-14.48609953*as.numeric(datMeth$cg13409216))+
  (-14.12054745*as.numeric(datMeth$cg13449372))+
  (0.148824189*as.numeric(datMeth$cg13460409))+
  (14.22126208*as.numeric(datMeth$cg13509147))+
  (1.146725133*as.numeric(datMeth$cg13510262))+
  (-13.34104874*as.numeric(datMeth$cg13514050))+
  (-16.79028681*as.numeric(datMeth$cg13550877))+
  (-7.65214571*as.numeric(datMeth$cg13564075))+
  (-2.123161461*as.numeric(datMeth$cg13571802))+
  (-0.096960459*as.numeric(datMeth$cg13587552))+
  (0.21020524*as.numeric(datMeth$cg13613532))+
  (23.83805501*as.numeric(datMeth$cg13631913))+
  (-0.764669972*as.numeric(datMeth$cg13654195))+
  (2.484283954*as.numeric(datMeth$cg13656062))+
  (14.36406636*as.numeric(datMeth$cg13656360))+
  (-5.203206279*as.numeric(datMeth$cg13700897))+
  (3.468680041*as.numeric(datMeth$cg13718960))+
  (0.043585555*as.numeric(datMeth$cg13843773))+
  (2.705537905*as.numeric(datMeth$cg13854874))+
  (0.588384291*as.numeric(datMeth$cg13861644))+
  (1.112549015*as.numeric(datMeth$cg13899108))+
  (3.583069518*as.numeric(datMeth$cg13975369))+
  (12.07884157*as.numeric(datMeth$cg13994175))+
  (1.041709284*as.numeric(datMeth$cg14009688))+
  (0.137366272*as.numeric(datMeth$cg14105047))+
  (-12.18111487*as.numeric(datMeth$cg14159818))+
  (-14.9517935*as.numeric(datMeth$cg14175438))+
  (-6.409663736*as.numeric(datMeth$cg14223995))+
  (21.11147672*as.numeric(datMeth$cg14281160))+
  (-4.233812041*as.numeric(datMeth$cg14350002))+
  (-4.795612884*as.numeric(datMeth$cg14423778))+
  (0.11617835*as.numeric(datMeth$cg14467840))+
  (-8.146118724*as.numeric(datMeth$cg14473016))+
  (-2.732393181*as.numeric(datMeth$cg14550518))+
  (-12.59980654*as.numeric(datMeth$cg14689355))+
  (7.024834167*as.numeric(datMeth$cg14747225))+
  (8.9044998*as.numeric(datMeth$cg14754581))+
  (7.981658933*as.numeric(datMeth$cg14916213))+
  (3.608249401*as.numeric(datMeth$cg14918082))+
  (0.385312405*as.numeric(datMeth$cg14972143))+
  (-2.533299173*as.numeric(datMeth$cg15013019))+
  (-1.981506525*as.numeric(datMeth$cg15171237))+
  (1.467052004*as.numeric(datMeth$cg15201877))+
  (-9.504209291*as.numeric(datMeth$cg15344028))+
  (0.002078841*as.numeric(datMeth$cg15381313))+
  (2.261768982*as.numeric(datMeth$cg15427448))+
  (17.06265852*as.numeric(datMeth$cg15447479))+
  (-0.925473409*as.numeric(datMeth$cg15489301))+
  (1.743616315*as.numeric(datMeth$cg15498283))+
  (5.216122264*as.numeric(datMeth$cg15551881))+
  (-12.95093372*as.numeric(datMeth$cg15569512))+
  (63.12415047*as.numeric(datMeth$cg15611364))+
  (0.032502686*as.numeric(datMeth$cg15642326))+
  (-8.009320965*as.numeric(datMeth$cg15811427))+
  (3.375380941*as.numeric(datMeth$cg15856055))+
  (-15.21848681*as.numeric(datMeth$cg15881088))+
  (-0.742899932*as.numeric(datMeth$cg15887846))+
  (-3.473480174*as.numeric(datMeth$cg15903282))+
  (-34.69842764*as.numeric(datMeth$cg15963417))+
  (-7.746394111*as.numeric(datMeth$cg15966757))+
  (-9.058118831*as.numeric(datMeth$cg16085042))+
  (8.01275733*as.numeric(datMeth$cg16173067))+
  (3.558270741*as.numeric(datMeth$cg16295988))+
  (6.023874287*as.numeric(datMeth$cg16313343))+
  (0.155277938*as.numeric(datMeth$cg16319578))+
  (32.14665692*as.numeric(datMeth$cg16340918))+
  (-4.990591369*as.numeric(datMeth$cg16354207))+
  (-1.880522626*as.numeric(datMeth$cg16357381))+
  (0.71058822*as.numeric(datMeth$cg16372520))+
  (4.024035102*as.numeric(datMeth$cg16408970))+
  (0.166685207*as.numeric(datMeth$cg16466334))+
  (-4.284958486*as.numeric(datMeth$cg16543027))+
  (2.080897193*as.numeric(datMeth$cg16612562))+
  (6.097618377*as.numeric(datMeth$cg16648841))+
  (3.187970744*as.numeric(datMeth$cg16713727))+
  (-4.410288204*as.numeric(datMeth$cg16718891))+
  (0.275798966*as.numeric(datMeth$cg16728114))+
  (3.265740439*as.numeric(datMeth$cg16743289))+
  (-0.761063169*as.numeric(datMeth$cg16816226))+
  (12.74522356*as.numeric(datMeth$cg16854606))+
  (7.416346316*as.numeric(datMeth$cg16933388))+
  (-1.359979716*as.numeric(datMeth$cg16984944))+
  (17.3157834*as.numeric(datMeth$cg17009433))+
  (-11.16561286*as.numeric(datMeth$cg17038116))+
  (9.679594816*as.numeric(datMeth$cg17129388))+
  (-8.786659598*as.numeric(datMeth$cg17133388))+
  (-4.043943288*as.numeric(datMeth$cg17324128))+
  (-0.762536285*as.numeric(datMeth$cg17431739))+
  (8.038276564*as.numeric(datMeth$cg17526300))+
  (-4.833656874*as.numeric(datMeth$cg17536848))+
  (-44.00939313*as.numeric(datMeth$cg17605084))+
  (1.069116234*as.numeric(datMeth$cg17627559))+
  (-9.926053105*as.numeric(datMeth$cg17641104))+
  (20.18800473*as.numeric(datMeth$cg17726022))+
  (6.468578944*as.numeric(datMeth$cg17749443))+
  (3.827665372*as.numeric(datMeth$cg17770886))+
  (0.023738294*as.numeric(datMeth$cg17861230))+
  (18.0100431*as.numeric(datMeth$cg17896249))+
  (0.058590013*as.numeric(datMeth$cg17903544))+
  (0.024808297*as.numeric(datMeth$cg17923358))+
  (-6.195165774*as.numeric(datMeth$cg17940013))+
  (-0.814304713*as.numeric(datMeth$cg17966192))+
  (-1.649720496*as.numeric(datMeth$cg18001427))+
  (-6.99391822*as.numeric(datMeth$cg18003795))+
  (11.53232389*as.numeric(datMeth$cg18117393))+
  (0.833493573*as.numeric(datMeth$cg18241647))+
  (13.32614877*as.numeric(datMeth$cg18267374))+
  (-0.960409715*as.numeric(datMeth$cg18384097))+
  (-1.196389134*as.numeric(datMeth$cg18392482))+
  (-0.514411149*as.numeric(datMeth$cg18468844))+
  (-0.137734976*as.numeric(datMeth$cg18587364))+
  (4.543967955*as.numeric(datMeth$cg18691434))+
  (-13.80959079*as.numeric(datMeth$cg18693704))+
  (0.797614179*as.numeric(datMeth$cg18732541))+
  (-0.430534396*as.numeric(datMeth$cg18771300))+
  (0.653578992*as.numeric(datMeth$cg18809289))+
  (-2.022128432*as.numeric(datMeth$cg18881501))+
  (33.05057141*as.numeric(datMeth$cg18996776))+
  (12.12001927*as.numeric(datMeth$cg19008809))+
  (-2.242872325*as.numeric(datMeth$cg19028160))+
  (-21.59866911*as.numeric(datMeth$cg19104072))+
  (-1.875110597*as.numeric(datMeth$cg19149785))+
  (-36.49384395*as.numeric(datMeth$cg19287114))+
  (-2.403650899*as.numeric(datMeth$cg19297232))+
  (-15.25553584*as.numeric(datMeth$cg19345165))+
  (-1.528350646*as.numeric(datMeth$cg19356189))+
  (0.402226134*as.numeric(datMeth$cg19371795))+
  (2.374343432*as.numeric(datMeth$cg19378133))+
  (35.83308293*as.numeric(datMeth$cg19398783))+
  (2.757897465*as.numeric(datMeth$cg19439331))+
  (3.775144564*as.numeric(datMeth$cg19514469))+
  (9.567877534*as.numeric(datMeth$cg19556572))+
  (2.560159507*as.numeric(datMeth$cg19560210))+
  (1.286269453*as.numeric(datMeth$cg19566405))+
  (-1.645326874*as.numeric(datMeth$cg19573166))+
  (5.285877086*as.numeric(datMeth$cg19586576))+
  (4.492646466*as.numeric(datMeth$cg19615059))+
  (1.89091587*as.numeric(datMeth$cg19632206))+
  (-2.827607692*as.numeric(datMeth$cg19663795))+
  (3.331314764*as.numeric(datMeth$cg19685066))+
  (-2.066714676*as.numeric(datMeth$cg19686152))+
  (-20.14577056*as.numeric(datMeth$cg19722847))+
  (-1.399107852*as.numeric(datMeth$cg19724470))+
  (2.914336066*as.numeric(datMeth$cg19731122))+
  (-3.218932622*as.numeric(datMeth$cg19883905))+
  (13.32700255*as.numeric(datMeth$cg20066677))+
  (4.874022045*as.numeric(datMeth$cg20090497))+
  (-12.17957599*as.numeric(datMeth$cg20162159))+
  (11.11358119*as.numeric(datMeth$cg20173259))+
  (1.829049562*as.numeric(datMeth$cg20234170))+
  (2.947623086*as.numeric(datMeth$cg20492933))+
  (-4.852301857*as.numeric(datMeth$cg20550118))+
  (-0.790590456*as.numeric(datMeth$cg20570279))+
  (3.062671547*as.numeric(datMeth$cg20572838))+
  (-2.37602848*as.numeric(datMeth$cg20652640))+
  (-0.120528788*as.numeric(datMeth$cg20674577))+
  (7.473261917*as.numeric(datMeth$cg20761322))+
  (2.193830696*as.numeric(datMeth$cg20828084))+
  (2.296781942*as.numeric(datMeth$cg20891917))+
  (5.85130181*as.numeric(datMeth$cg20967028))+
  (3.152813941*as.numeric(datMeth$cg21006686))+
  (-0.052460218*as.numeric(datMeth$cg21053529))+
  (-2.970173829*as.numeric(datMeth$cg21081971))+
  (-7.319261619*as.numeric(datMeth$cg21099326))+
  (-5.286945909*as.numeric(datMeth$cg21120249))+
  (-5.983252062*as.numeric(datMeth$cg21137706))+
  (-1.102666731*as.numeric(datMeth$cg21184495))+
  (3.959178135*as.numeric(datMeth$cg21200703))+
  (0.648066076*as.numeric(datMeth$cg21201109))+
  (6.268763011*as.numeric(datMeth$cg21207418))+
  (15.06457901*as.numeric(datMeth$cg21296230))+
  (1.140495095*as.numeric(datMeth$cg21363706))+
  (-1.462031413*as.numeric(datMeth$cg21649520))+
  (-1.765357537*as.numeric(datMeth$cg21712685))+
  (2.116884479*as.numeric(datMeth$cg21762589))+
  (19.61967861*as.numeric(datMeth$cg21801378))+
  (3.158955159*as.numeric(datMeth$cg21835643))+
  (2.497611714*as.numeric(datMeth$cg21907579))+
  (-2.151586582*as.numeric(datMeth$cg21926612))+
  (1.01756345*as.numeric(datMeth$cg21993406))+
  (-0.38718029*as.numeric(datMeth$cg22090592))+
  (-1.509071514*as.numeric(datMeth$cg22179082))+
  (0.727265735*as.numeric(datMeth$cg22194129))+
  (-0.117965431*as.numeric(datMeth$cg22197830))+
  (4.813692319*as.numeric(datMeth$cg22282672))+
  (0.027494572*as.numeric(datMeth$cg22395019))+
  (6.4288385*as.numeric(datMeth$cg22396353))+
  (4.232471361*as.numeric(datMeth$cg22407458))+
  (4.044042377*as.numeric(datMeth$cg22473095))+
  (1.547278123*as.numeric(datMeth$cg22484793))+
  (9.772363822*as.numeric(datMeth$cg22495124))+
  (5.567558371*as.numeric(datMeth$cg22511262))+
  (6.279600362*as.numeric(datMeth$cg22512531))+
  (3.485338581*as.numeric(datMeth$cg22580353))+
  (6.352315*as.numeric(datMeth$cg22582569))+
  (16.6352922*as.numeric(datMeth$cg22594309))+
  (31.84219355*as.numeric(datMeth$cg22736354))+
  (1.377506045*as.numeric(datMeth$cg22809047))+
  (-7.872570922*as.numeric(datMeth$cg22947000))+
  (8.613889336*as.numeric(datMeth$cg22971191))+
  (0.434401518*as.numeric(datMeth$cg22983092))+
  (-5.681543236*as.numeric(datMeth$cg22991148))+
  (-5.804862482*as.numeric(datMeth$cg23124451))+
  (1.686038152*as.numeric(datMeth$cg23127998))+
  (2.580749951*as.numeric(datMeth$cg23152772))+
  (-0.632208093*as.numeric(datMeth$cg23159337))+
  (-0.047747114*as.numeric(datMeth$cg23173910))+
  (3.958521328*as.numeric(datMeth$cg23191950))+
  (1.448091683*as.numeric(datMeth$cg23213217))+
  (2.789697779*as.numeric(datMeth$cg23234999))+
  (0.774552515*as.numeric(datMeth$cg23239039))+
  (1.264378979*as.numeric(datMeth$cg23338195))+
  (-7.470945638*as.numeric(datMeth$cg23376526))+
  (-0.644984918*as.numeric(datMeth$cg23568913))+
  (-3.15460847*as.numeric(datMeth$cg23668631))+
  (24.15870289*as.numeric(datMeth$cg23710218))+
  (3.603684682*as.numeric(datMeth$cg23818978))+
  (-31.90246973*as.numeric(datMeth$cg23832061))+
  (-1.01239731*as.numeric(datMeth$cg24110063))+
  (7.120190402*as.numeric(datMeth$cg24125648))+
  (22.54068526*as.numeric(datMeth$cg24208206))+
  (23.3003573*as.numeric(datMeth$cg24304712))+
  (0.617734017*as.numeric(datMeth$cg24332433))+
  (-12.95513488*as.numeric(datMeth$cg24407308))+
  (-4.202135081*as.numeric(datMeth$cg24493940))+
  (0.307545432*as.numeric(datMeth$cg24505122))+
  (2.841188478*as.numeric(datMeth$cg24505341))+
  (1.410756917*as.numeric(datMeth$cg24556026))+
  (0.653433186*as.numeric(datMeth$cg24651706))+
  (5.41797411*as.numeric(datMeth$cg24674703))+
  (-4.986411499*as.numeric(datMeth$cg24921089))+
  (-0.085880587*as.numeric(datMeth$cg25022327))+
  (10.80113727*as.numeric(datMeth$cg25092328))+
  (-9.522144355*as.numeric(datMeth$cg25136687))+
  (5.410729068*as.numeric(datMeth$cg25229964))+
  (2.248622897*as.numeric(datMeth$cg25251635))+
  (-9.36844987*as.numeric(datMeth$cg25256723))+
  (-13.6572112*as.numeric(datMeth$cg25428451))+
  (-5.900533275*as.numeric(datMeth$cg25459323))+
  (-2.837642113*as.numeric(datMeth$cg25536676))+
  (2.25600198*as.numeric(datMeth$cg25713185))+
  (-7.808838536*as.numeric(datMeth$cg25769980))+
  (-15.99899499*as.numeric(datMeth$cg25881193))+
  (-0.725019577*as.numeric(datMeth$cg25898500))+
  (-2.270877875*as.numeric(datMeth$cg26022315))+
  (-17.94903039*as.numeric(datMeth$cg26091688))+
  (7.269147318*as.numeric(datMeth$cg26096837))+
  (-5.768101285*as.numeric(datMeth$cg26104204))+
  (-3.984563973*as.numeric(datMeth$cg26109803))+
  (5.446299477*as.numeric(datMeth$cg26201213))+
  (-8.269454924*as.numeric(datMeth$cg26212924))+
  (1.046570867*as.numeric(datMeth$cg26219051))+
  (-1.182810193*as.numeric(datMeth$cg26312920))+
  (-1.261236982*as.numeric(datMeth$cg26350286))+
  (-10.00204065*as.numeric(datMeth$cg26357744))+
  (40.42085373*as.numeric(datMeth$cg26382071))+
  (7.25927706*as.numeric(datMeth$cg26394737))+
  (-0.852196304*as.numeric(datMeth$cg26394940))+
  (-3.998837847*as.numeric(datMeth$cg26581729))+
  (-2.026512819*as.numeric(datMeth$cg26614073))+
  (-11.15430006*as.numeric(datMeth$cg26665419))+
  (2.558421247*as.numeric(datMeth$cg26711820))+
  (3.866139002*as.numeric(datMeth$cg26746469))+
  (-17.14738745*as.numeric(datMeth$cg26815229))+
  (4.586931755*as.numeric(datMeth$cg26824091))+
  (0.063834719*as.numeric(datMeth$cg26842024))+
  (9.372871086*as.numeric(datMeth$cg26866325))+
  (-2.24450608*as.numeric(datMeth$cg26898166))+
  (-3.259148427*as.numeric(datMeth$cg26932976))+
  (-16.80473407*as.numeric(datMeth$cg27015931))+
  (-33.54555855*as.numeric(datMeth$cg27187881))+
  (-0.031981428*as.numeric(datMeth$cg27244482))+
  (-15.82375956*as.numeric(datMeth$cg27367952))+
  (0.033980159*as.numeric(datMeth$cg27440834))+
  (-21.2007229*as.numeric(datMeth$cg27493997))+
  (-4.390952425*as.numeric(datMeth$cg27514224))+
  (-1.465742937*as.numeric(datMeth$cg27626102))+
  (9.950424372*as.numeric(datMeth$cg27655905))


####Extract these 500 cpgs from my matrix, then if there is a missing value in the column then replace with a mean of the methylation value (ie the row).
##with calculating the mean, use na.rm, eg mean(cpg1, na.rm =T)

##adding in age to meta
GSE61461.metaDat2$DNAmPhenoAge <- DNAmPhenoAge



###############Weidner Wolfgang 99 CpG Clock  ####################################


datMeth <- t(read.csv("GSE614161_Fullarray.csv", header = T, row.names = 1))

WW99 <- read.csv("../../Project Rejuvination/R_dir_P.Rejuv/weidner_99_CpGs_coef.csv", header=T)


cpgdat <- datMeth[,as.character(WW99$CpG.sites[-1])]


##list of coef without the intercept value
conew <- WW99$Coefficient[-1]

#Intercept standard deviation, the first big value that all the other values are subracted from.
sumest <- WW99$Coefficient[1]

#####Imputing NAs as zeros in new object, you dont have to do this, this was something I happened to do here
#GSE54848.99.na0 <- GSE54848.99[,as.character(WW99$CpG.sites)[-1]]

#GSE54848.99.na0[is.na(GSE54848.99.na0)] <- 0

#cpgdat <- GSE54848.99.na0

#########Age predictor function ##################################################
predinator<-function(demdats){
  predAge=0
  newval=0
  count=0
  for(val in demdats){
    #print(val)
    count=count+1
    num<- (conew[count]*val)
    newval<-newval+num
    #print(count)
  }
  #adding or subtracting values from this line can correct age prediction sometimes
  predAge<-sumest + newval 
  print(predAge)
}

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

Weidnerpred99 <- as.numeric(predout)


##adding in age to meta
GSE61461.metaDat2$Weidnerpred99 <- Weidnerpred99

############## 3 CpG clock ##################

WW3 <- read.csv("../../Project Rejuvination/R_dir_P.Rejuv/3_CpG_Weidner.csv", header=T)


cpgdat <- datMeth[,as.character(WW3$CpG.sites[-1])]

###just using days as actual age here
#ages <- days2

##list of coef without the intercept value
conew <- WW3$Coefficient[-1]

#Intercept standard deviation, the first big value that all the other values are subracted from.
sumest <- WW3$Coefficient[1]



#########Age predictor function ##################################################
predinator<-function(demdats){
  predAge=0
  newval=0
  count=0
  for(val in demdats){
    #print(val)
    count=count+1
    num<- (conew[count]*val)
    newval<-newval+num
    #print(count)
  }
  #adding or subtracting values from this line can correct age prediction sometimes
  predAge<-sumest + newval 
  print(predAge)
}

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

Weidnerpred3 <- as.numeric(predout)

##adding in age to meta
GSE61461.metaDat2$Weidnerpred3 <- Weidnerpred3


#############Hannum age #######################


hannum = read.table("../../Project Rejuvination/R_dir_P.Rejuv/Hannum_etal_age_predictive_probes.txt", header=T, sep="\t")
hannum <- hannum[,c(1,6)]
rownames(hannum) <- hannum$Marker
cpgdat <- t(datMeth)
probes <- intersect(hannum$Marker, rownames(cpgdat))

b = cpgdat[probes,]
p = hannum[probes,]

for (i in probes) {
  b[i,]= as.numeric(b[i,])*as.numeric(p[i,"Coefficient"])
}

predicted_age=colSums(b)
hannum_pred_age=predicted_age
hannum_pred_age <- as.data.frame(hannum_pred_age)


##adding in age to meta
GSE61461.metaDat2$hannum_pred_age <- hannum_pred_age$hannum_pred_age


###Teloage #########

Telo <- read.csv(file.choose(), header=T)


cpgdat <- datMeth[,as.character(Telo$Name[-1])]

###this is so weird, only one CpG present...doesnt make sense.

###just using days as actual age here
#ages <- days2

##list of coef without the intercept value
conew <- Telo$Coeficient[-1]

#Intercept standard deviation, the first big value that all the other values are subracted from.
sumest <- Telo$Coeficient[1]



#########Age predictor function ##################################################
predinator<-function(demdats){
  predAge=0
  newval=0
  count=0
  for(val in demdats){
    #print(val)
    count=count+1
    num<- (conew[count]*val)
    newval<-newval+num
    #print(count)
  }
  #adding or subtracting values from this line can correct age prediction sometimes
  predAge<-sumest + newval 
  print(predAge)
}

####for predicting 
predout=0
for(x in 1:nrow(cpgdat)){
  #print(x)
  predout[x] <- predinator(cpgdat[x,])
}

TeloPred <- predout

##adding in age to meta
GSE61461.metaDat2$TeloPred <- TeloPred


####Example plot ######

#horvath age

#temp dataframe
big1 <- data.frame(HorvathAge = GSE61461.metaDat2$Horvath,
                   ActualAge = GSE61461.metaDat2$age,
                   Status = GSE61461.metaDat2$genotype)

#For other plots, just change Horvath to any other clocks youve added to your metadata object. Change status to anything else, say father/son? Remove if not needed.

RdYlBu <- brewer.pal(n = 9, name = "RdYlBu")

cols <- c("darkgrey", "#FDAE61","#ABD9E9")                   
pdf("Horvath-lm.pdf", width=7, height=5, useDingbats=FALSE)
ggplot(big1, aes(x=ActualAge, y=HorvathAge, color = Status))  + labs(y = "Horvath Age", x = "Chronological Age") +
  scale_color_manual(values = cols) + theme_classic() + geom_point()  + geom_smooth(method="lm", fill = NA)
dev.off() 



