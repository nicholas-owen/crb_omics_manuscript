library("methylKit")

directory <- "../01_Preprocess/NUDUP/Ready/"

files <- dir(directory)
# This line is only used to extract Sample Name.
file.name <- as.list(unname(sapply(files,function(x) strsplit(x,split="_")[[1]][1])))
files <- as.list(paste(directory,files,sep=""))

message("Processing Bismark Alignment...")
my.methRaw=processBismarkAln(location=files, 
                             sample.id=file.name, 
                             assembly="danRer11", 
                             save.context=NULL, 
                             read.context="CpG", 
                             mincov=7,
                             minqual=20,
                             treatment=c(0,0,0,1,1,1),
                             save.folder="./")

message("Reading CpG into R with methylKit...")
if (!file.exists("./Data")) dir.create("./Data")
CpGFiles <- as.list(paste("./",unlist(file.name),"_CpG.txt",sep=""))
sampleName <- c("WildType_1", "WildType_2", "WildType_3", "Mutation_1", "Mutation_2", "Mutation_3")

myobj <- methRead(CpGFiles,
                  sample.id= as.list(sampleName),
                  assembly="danRer11",
                  treatment=c(0,0,0,1,1,1),
                  context="CpG",
                  mincov=7)

save(myobj, file="./Data/myobj.rda")

message("Figure: plot distribution")
if (!file.exists("./Figure")) dir.create("./Figure")

pdf("./Figure/DensityPlot.pdf", width=16, height=8)
par(mfrow=c(2,3))
for(i in 1:length(myobj)) getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
dev.off()


message("Figure: plot coverage")
if (!file.exists("./Figure")) dir.create("./Figure")

pdf("./Figure/CoveragePlot.pdf", width=16, height=8)
par(mfrow=c(2,3))
for(i in 1:length(myobj)) getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
dev.off()


message("Filtering Coverage and Normalisation")
filtered.myobj <- filterByCoverage(myobj,lo.count=7,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
filtered.normed.myobj <- normalizeCoverage(filtered.myobj)

message("Figure: plot filtered coverage")
if (!file.exists("./Figure")) dir.create("./Figure")

pdf("./Figure/FilteredCoveragePlot.pdf", width=16, height=8)
par(mfrow=c(2,3))
for(i in 1:length(filtered.normed.myobj)) getCoverageStats(filtered.normed.myobj[[i]], plot=TRUE, both.strands=FALSE)
dev.off()

message("Create Meth and Beta Object")
if (!file.exists("./Data")) dir.create("./Data")

meth <- unite(filtered.normed.myobj, destrand=FALSE)
beta <- percMethylation(meth)/100

save(meth,file="./Data/meth.rda")
save(beta,file="./Data/beta.rda")

