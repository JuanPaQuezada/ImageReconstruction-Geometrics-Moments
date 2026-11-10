library(jpeg)
library(ggplot2)
library(EBImage)
library(gridExtra)
img1<-readJPEG("/home/juanpa333/Documents/ImageReconstruction-Geometrics-Moments/data/raw/Experimento4_par1.jpg")
#img1<-readJPEG(file.choose())
matrix_gris<-function(image){
  #convert image to grayscale if it is not already
  if (length(dim(image)) == 3) {
    return ((image[,,1]+image[,,2]+image[,,3])/3)
  }
  return(image)
}

imagen_gris<-function(matrix_gris){
  #display the grayscale image
  display(matrix_gris, method = "raster")
  return(matrix_gris)
}

otsu_threshold<-function(image){
    
    matrix_gray=matrix_gris(image)
    vector_gray=as.numeric(matrix_gray)
    if(max(vector_gray)<=1){
        vector_gray=vector_gray*255
    }
    #calculate histogram of the image 
    histo<-hist(vector_gray, breaks=0:256, plot=TRUE, main="Histograma de la imagen", xlab="Niveles de gris", ylab="Frecuencia") 
    gray_levels<-histo$mids
    counts<-histo$counts
    imagen_raster<-as.raster(matrix_gray)
    plot(imagen_raster,main="Imagen")
    grey_levels<-histo$mids
    counts<-histo$counts
    total_pixels<-sum(counts)
    prob<-counts/total_pixels
    #histogram of levels of gray 
    hist(prob, breaks=seq(0,1,by=0.01), main="Histograma de niveles de gris", xlab="Niveles de gris", ylab="Frecuencia")
    #calculate cumulative sum of probabilities and cumulative mean 
    omega<-cumsum(prob)
    mu<-cumsum(prob*gray_levels)
    mu_t<-mu[length((mu))]
    sigma_b_squared<-(mu_t*omega-mu)^2/(omega*(1-omega))
    sigma_b_squared[is.nan(sigma_b_squared)]<-0
    threshold_index<-which.max(sigma_b_squared)
    threshold_value<-gray_levels[threshold_index]
    cat("El umbral de Otsu es:", threshold_value, "\n")
    umbral_normalized<-threshold_value/255
    cat("El umbral de Otsu normalizado es:", umbral_normalized, "\n")
    return(umbral_normalized)
}
matriz_gris<-matrix_gris(img1)
img_gray<-imagen_gris(matriz_gris)
otsu_threshold(img_gray)
#mostrar imagen en escala de grises

