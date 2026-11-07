
# ============================================================
# REGISTRO DE IMAGENES BINARIAS MEDIANTE MOMENTOS
#
# 1. Lectura de imagenes
# 2. Umbralizacion
# 3. Centrado (traslacion)
# 4. Orientacion (rotacion)
# 5. Ajuste de tamano (escalamiento)
#
# La imagen 1 es la referencia.
# Solamente se transforma la imagen 2.
# ============================================================

# ------------------------------------------------------------
# 0. PAQUETES
# ------------------------------------------------------------

library(jpeg)

# ------------------------------------------------------------
# 1. LEER LAS IMAGENES
# ------------------------------------------------------------

img1 <- readJPEG(file.choose())
img2 <- readJPEG(file.choose())

# Convertir a escala de grises
gris <- function(img) {
  
  if (length(dim(img)) == 3) {
    return( (img[,,1] + img[,,2] + img[,,3])/3 )
  }
  
  img
}

g1 <- gris(img1)
g2 <- gris(img2)

# ------------------------------------------------------------
# 2. UMBRALIZACION
# ------------------------------------------------------------

# El objeto es mas oscuro que el fondo

umbral1 <- 0.4  # Estos valores se deben ajustar para cada imagen
umbral2 <- 0.4  # Estos valores se deben ajustar para cada imagen

b1 <- g1 < umbral1 # Esto es para objeto obscuro, fondo claro
b2 <- g2 < umbral2 # Si es objeto claro, fondo obscuro, se cambian a >

# TRUE = objeto
# FALSE = fondo

# ------------------------------------------------------------
# 3. FUNCION PARA CALCULAR MOMENTOS
# ------------------------------------------------------------

momentos <- function(B) {
  
  pix <- which(B, arr.ind = TRUE)
  
  if (nrow(pix) == 0)
    stop("La imagen no contiene pixeles del objeto.")
  
  # Coordenadas
  y <- pix[,1]
  x <- pix[,2]
  
  # Momento de orden cero
  M00 <- length(x)
  
  # Centroide
  xc <- sum(x) / M00
  yc <- sum(y) / M00
  
  # Momentos centrales de segundo orden
  mu20 <- sum((x - xc)^2)
  mu02 <- sum((y - yc)^2)
  mu11 <- sum((x - xc)*(y - yc))
  
  # Orientacion principal
  theta <- 0.5 * atan2(
    2 * mu11,
    mu20 - mu02
  )
  
  list(
    area = M00,
    xc = xc,
    yc = yc,
    mu20 = mu20,
    mu02 = mu02,
    mu11 = mu11,
    theta = theta
  )
}

# ------------------------------------------------------------
# 4. FUNCION PARA TRANSFORMAR UNA IMAGEN
# ------------------------------------------------------------

# Utiliza mapeo inverso y vecino mas cercano.
#
# centro_mov: centroide en la imagen original
# centro_ref: centroide deseado en la imagen resultante
# angulo: rotacion en radianes
# escala: factor de escalamiento

transformar <- function(B, ref,
                        centro_mov,
                        centro_ref,
                        angulo = 0,
                        escala = 1) {
  
  nr <- nrow(ref)
  nc <- ncol(ref)
  
  # Coordenadas de la imagen resultante
  coords <- expand.grid(
    x = 1:nc,
    y = 1:nr
  )
  
  # Coordenadas relativas al centro de referencia
  dx <- coords$x - centro_ref[1]
  dy <- coords$y - centro_ref[2]
  
  # Transformacion inversa
  x0 <- (cos(angulo)*dx + sin(angulo)*dy) / escala
  
  y0 <- (-sin(angulo)*dx + cos(angulo)*dy) / escala
  
  # Regresar al sistema original
  x0 <- x0 + centro_mov[1]
  y0 <- y0 + centro_mov[2]
  
  # Vecino mas cercano
  xi <- round(x0)
  yi <- round(y0)
  
  # Imagen resultante
  salida <- matrix(
    FALSE,
    nrow = nr,
    ncol = nc
  )
  
  # Coordenadas validas
  validos <- xi >= 1 & xi <= ncol(B) &
    yi >= 1 & yi <= nrow(B)
  
  # Asignar los pixeles
  salida[cbind(
    coords$y[validos],
    coords$x[validos]
  )] <- B[cbind(
    yi[validos],
    xi[validos]
  )]
  
  salida
}

# ------------------------------------------------------------
# 5. CALCULAR MOMENTOS INICIALES
# ------------------------------------------------------------

m1 <- momentos(b1)
m2 <- momentos(b2)

cat("\nMOMENTOS INICIALES\n")

print(m1)
print(m2)

# ============================================================
# 6. PRIMERA TRANSFORMACION: TRASLACION
# ============================================================

# Calcular desplazamiento necesario
tx <- m1$xc - m2$xc
ty <- m1$yc - m2$yc

cat("\nTRASLACION\n")
cat("tx =", tx, "\n")
cat("ty =", ty, "\n")

# Aplicar traslacion
b2_tras <- transformar(
  B = b2,
  ref = b1,
  centro_mov = c(m2$xc, m2$yc),
  centro_ref = c(m1$xc, m1$yc)
)

# ------------------------------------------------------------
# RECALCULAR MOMENTOS DESPUES DE TRASLADAR
# ------------------------------------------------------------

m2_tras <- momentos(b2_tras)

cat("\nMOMENTOS DESPUES DE TRASLADAR\n")
print(m2_tras)

# ============================================================
# 7. SEGUNDA TRANSFORMACION: ROTACION
# ============================================================

# Diferencia entre orientaciones
dtheta <- m1$theta - m2_tras$theta

# Resolver ambiguedad de 180 grados
dtheta <- (dtheta + pi/2) %% pi -pi/2   # a veces es  +pi/2

cat("\nROTACION\n")
cat("Angulo (radianes) =", dtheta, "\n")
cat("Angulo (grados) =", dtheta*180/pi, "\n")

# Aplicar rotacion
# Se utiliza el centroide actualizado de la imagen trasladada
# y el centroide de referencia como destino.

b2_rot <- transformar(
  B = b2_tras,
  ref = b1,
  centro_mov = c(m2_tras$xc, m2_tras$yc),
  centro_ref = c(m1$xc, m1$yc),
  angulo = dtheta
)

# ------------------------------------------------------------
# RECALCULAR MOMENTOS DESPUES DE ROTAR
# ------------------------------------------------------------

m2_rot <- momentos(b2_rot)

cat("\nMOMENTOS DESPUES DE ROTAR\n")
print(m2_rot)

# ============================================================
# 8. TERCERA TRANSFORMACION: ESCALAMIENTO
# ============================================================

# Calcular factor de escala utilizando las areas actualizadas
s <- sqrt(m1$area / m2_rot$area)

cat("\nESCALAMIENTO\n")
cat("Factor de escala =", s, "\n")

# Aplicar escalamiento
# Se centra en el centroide actualizado y se lleva
# al centroide de referencia.

b2_reg <- transformar(
  B = b2_rot,
  ref = b1,
  centro_mov = c(m2_rot$xc, m2_rot$yc),
  centro_ref = c(m1$xc, m1$yc),
  escala = s
)

# ------------------------------------------------------------
# RECALCULAR MOMENTOS FINALES
# ------------------------------------------------------------

m2_final <- momentos(b2_reg)

cat("\nMOMENTOS FINALES\n")
print(m2_final)

# ============================================================
# 9. RESUMEN DE RESULTADOS
# ============================================================

cat("\n====================================\n")
cat("RESUMEN DEL REGISTRO\n")
cat("====================================\n")

cat("\nTraslacion:\n")
cat("tx =", tx, "\n")
cat("ty =", ty, "\n")

cat("\nRotacion:\n")
cat("Grados =", dtheta*180/pi, "\n")

cat("\nEscalamiento:\n")
cat("s =", s, "\n")

cat("\nAreas:\n")
cat("Referencia =", m1$area, "\n")
cat("Registrada =", m2_final$area, "\n")

cat("\nCentroides:\n")
cat("Referencia =", m1$xc, m1$yc, "\n")
cat("Registrada =", m2_final$xc, m2_final$yc, "\n")

# ============================================================
# 10. VISUALIZACION DE LAS TRANSFORMACIONES
# ============================================================

mostrar <- function(B, titulo) {
  
  image(
    t(B[nrow(B):1, ]),
    col = c("white", "black"),
    axes = FALSE,
    main = titulo
  )
}

par(mfrow = c(2,3), mar = c(2,2,3,1))

mostrar(b1, "Imagen 1: referencia")
mostrar(b2, "Imagen 2: original")
mostrar(b2_tras, "Traslacion")
mostrar(b2_rot, "Rotacion")
mostrar(b2_reg, "Escalamiento")

# ------------------------------------------------------------
# SUPERPOSICION FINAL
# ------------------------------------------------------------

# Rojo: imagen de referencia
# Azul: imagen registrada
# Magenta: coincidencia

superposicion <- array(
  0,
  dim = c(nrow(b1), ncol(b1), 3)
)

superposicion[,,1] <- b1
superposicion[,,3] <- b2_reg

plot(
  0, 0,
  type = "n",
  xlim = c(0,1),
  ylim = c(-1,0),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Superposicion final"
)

rasterImage(
  as.raster(superposicion[nrow(b1):1,,]),
  0, 0, 1, -1
)

par(mfrow = c(1,1))

# ============================================================
# Visualización simplificada para el reporte
# ============================================================
par(mfrow = c(2,2), mar = c(2,2,3,1))

mostrar(b1, "Imagen 1: referencia")
mostrar(b2, "Imagen 2: original")
mostrar(b2_reg, "Registo")

# ------------------------------------------------------------
# SUPERPOSICION FINAL
# ------------------------------------------------------------

# Rojo: imagen de referencia
# Azul: imagen registrada
# Magenta: coincidencia

superposicion <- array(
  0,
  dim = c(nrow(b1), ncol(b1), 3)
)

superposicion[,,1] <- b1
superposicion[,,3] <- b2_reg

plot(
  0, 0,
  type = "n",
  xlim = c(0,1),
  ylim = c(-1,0),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Superposicion final"
)

rasterImage(
  as.raster(superposicion[nrow(b1):1,,]),
  0, 0, 1, -1
)

par(mfrow = c(1,1))

# ============================================================
# ============================================================

