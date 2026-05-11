#author: Pedro Carrillo Alarcon
#Carga de librerias
library (org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(statmod)

#Set work directory
setwd(dir = "/home/pcarrilloalarcon/transcriptomica/Proyecto_2026/")
#load raw count matrix
seqdata = read.delim(file = "Results/GSE60450_Lactation-GenewiseCounts.txt",
                     row.names = "EntrezGeneID")
seqdata = seqdata[,2:ncol(seqdata)]

# low Metadatos y determinar groups
sampleinfo = read.delim(file = "metadata.txt")
group = paste(sampleinfo$CellType, sampleinfo$Status, sep = ".")
group = as.factor(group)


#Rename del dataset
colnames(seqdata) = substr(x = colnames(seqdata), start = 1, stop = 7)

#Obtener valores maximos por gen del dataset y los genes asocaidos a estos valores
valores_max = sapply(seqdata, max, na.rm = TRUE)
gen_asociado <- apply(seqdata, 2, function(x) rownames(seqdata)[which.max(x)] )

#Crea un dataframe con el valor de los nombres de las muestras, el valor maximo y el gen asociado a estos valores
muestras = colnames(seqdata)
df = data.frame(
  muestra = muestras,
  valores_maximos = valores_max,
  gen_asociado = gen_asociado
)
#Greáfico expresando el df anterior
ggplot(df, mapping= aes(x=muestra, y=valores_maximos)) + geom_col()+ labs(title = "Valores Máximos por Muestra")

#Calcular el porcentaje de genes que tiene un valor > 0 por muestra
porcentaje_genes <- apply(seqdata > 0, 2, sum) / nrow(seqdata) * 100

#Añadir el porcentaje de genes al df
df = data.frame(
  muestra = muestras,
  valores_maximos = valores_max,
  gen_asociado = gen_asociado,
  pocentaje_genes_mayor_0 = porcentaje_genes
)

#Transformar los IDs con el ID oriinal, el simbolo y el nombre completo
ann <- select(org.Mm.eg.db,keys=rownames(seqdata),
              columns=c("ENTREZID","SYMBOL","GENENAME"))
head(ann)
#Determianmos cuales de los GENID no han podido actualizar al nuevo formato
na_values = sum(is.na(ann$SYMBOL))

#TRansformar a objeto DGE
y <- DGEList(seqdata)

#Determinar el valor de las librerias (Deben ser tamaños grandes, si el exp va bien)
sum(colSums(seqdata))

#Añade al objeto DGE los grupos determinados de los metadatos
y$samples$group = group
y$genes = ann

#Filtrar por low-expression
keep = filterByExpr(y)
summary(keep)

#False para que vuelva a calcular el tamaño de la libreria, eliminando los falsos del filtrado
y <- y[keep, , keep.lib.sizes=FALSE]
head (y)

#Transformar los datos para que sean legibles
cpm = cpm(y, log=TRUE)

#Grafico pre-normalizacion
boxplot(cpm, las=2, main="LogCPM por muestra (sin normalizar)", 
        ylab="logCPM", cex.axis=0.8)

#Normalizacion de los datos
y = calcNormFactors(y)

#Aislar el valor de los datos normalizados
norm_factors = y$samples$norm.factors

#Determinar que muestras tiene un valor de normalizacion <1
muestras_bajo1 = rownames(y$samples)[norm_factors < 1]

#Mostrar las muestras con este nivel
print(muestras_bajo1)

#Volver a transfroamr los datos tras la normlaizacion
cpm = cpm(y, log=TRUE)
#Gráfico post-normalizacion
boxplot(cpm, las=2, main="LogCPM por muestra (normalizado)", 
        ylab="logCPM", cex.axis=0.8)

#Generar gráfico MDS para las muestras
pch <- c(0,1,2,15,16,17)

colors <- rep(c("darkgreen", "red", "blue"), 2)

plotMDS(y, col=colors[group], pch=pch[group])

legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2, cex = 0.8)

#Crear la matriz de diseño
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

#Calculo de valores de tipo de dispersion
val_dispersion = estimateDisp(y, design, robust=TRUE)

#Crea vector con los nombres de los genes que quieres calcular
nombres = c("Imp4", "Xkr4", "Sox17")
# Esto te devolverá  el valor de tagwise de Xkr4 Imp4 y Sox17 y lo guardará en un vector
tagwise = c(val_dispersion$tagwise.dispersion[which(val_dispersion$genes$SYMBOL == "Imp4")],
val_dispersion$tagwise.dispersion[which(val_dispersion$genes$SYMBOL == "Xkr4")],
val_dispersion$tagwise.dispersion[which(val_dispersion$genes$SYMBOL == "Sox17")])

# Esto te devolverá  el valor de tagwise de Xkr4 Imp4 y Sox17 y lo guardará en un vector
trended = c(val_dispersion$trended.dispersion[which(val_dispersion$genes$SYMBOL == "Imp4")],
            val_dispersion$trended.dispersion[which(val_dispersion$genes$SYMBOL == "Xkr4")],
            val_dispersion$trended.dispersion[which(val_dispersion$genes$SYMBOL == "Sox17")])
#Crea un df para almacenar todos los datos de forma ordenada
df_dispersion = data.frame(genes = nombres,
                          tagwise_dispersion = tagwise,
                          trended_dispersion = trended)

modelo_GLM = glmQLFit(val_dispersion, design, robust=TRUE)

#Calcular coeficientes para genes mencionados anteriormente
#Guarda el valor de la fila que equivale a cada gen para buscarlo en el objeto DGE
indices <- which(modelo_GLM$genes$SYMBOL %in% nombres)
#Guarda el valor de los coeficientes asociados a cada valor indice
coeficientes_GLM = modelo_GLM$coefficients[indices, ]
#Cambia el nombre del objeto creado anteriormente para que sean los nombres de los genes los que aparezcan como rownames
rownames(coeficientes_GLM) = y$genes$SYMBOL[indices]

#Prueba de hipótesis
#definir la comparacion
B.LvsP <- makeContrasts(basal.lactate-basal.pregnant, levels=design)

#Realiza la prueba de hipótesis con la condicion anteriormente definida, obteniendo logFc o Pvalue
res <- glmQLFTest(modelo_GLM, contrast= B.LvsP )
#head(res$table)
#Ordena los reusltados y calcula el FDR
resultados = topTags(res, n=Inf)
#dim(resultados)
#head(resultados$table)

#Filtrar por el valor de FDR, indica los genes que cumplen con la condicion.
nrow(resultados$table[resultados$table$FDR <= 0.05, ])
nrow(resultados$table[resultados$table$PValue <= 0.05, ])


#Filtra los genes de la prueba topTags que tiennen un pvalue menor a 0.05 y un fc mayor a uno. R ya sabe que es mayor o gual o menor o igual, por eso solo se pone igual.
res_degs <- decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 1)
head(res_degs)
#Realiza un conteo rápido y organizado de cuantos genes han pasado los filtros de significancia. NotSig, son los que no y el resto si.
summary (res_degs)




