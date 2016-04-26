#### Divide Testing and Training Data

#Stripes
test_raster_A <- plot_boundary
values(test_raster_A) <- rep(0:1, each = 40, length.out = 402*166)
plot(test_raster_A)

# Better Stripes --- DO THIS ONE!
test_raster_AA <- plot_boundary
values(test_raster_AA) <- c(rep(c(1,1,0,1,0,1), each = 67, length.out = 402*41),
                            rep(c(0,1,1,1,0,1), each = 67, length.out = 402*42),
                            rep(c(0,1,0,1,1,1), each = 67, length.out = 402*41),
                            rep(c(0,1,0,1,0,1), each = 67, length.out = 402*42))
plot(test_raster_AA)

# 0 = test
# 1 = train


#### Remove 1/3
test_raster_B <- plot_boundary
values(test_raster_B) <- rep(c(0,1,1), each = 134, length.out = 402*166)
plot(test_raster_B)

#### Checkerboard - EVEN
test_raster_C <- plot_boundary
values(test_raster_C) <- c(rep(c(0,1,0,0,1,0), each = 134, length.out = 402*83),
                           rep(c(1,0,1,1,0,1), each = 134, length.out = 402*83))
plot(test_raster_C)

#### Checkerboard - UN-EVEN
test_raster_D <- plot_boundary
values(test_raster_D) <- c(rep(c(0,1,0,1,0,1), each = 67, length.out = 402*41),
                           rep(c(1,0,1,0,1,0), each = 67, length.out = 402*42),
                           rep(c(0,1,0,1,0,1), each = 67, length.out = 402*41),
                           rep(c(1,0,1,0,1,0), each = 67, length.out = 402*42))
plot(test_raster_D)



