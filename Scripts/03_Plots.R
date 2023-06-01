##### [F1] Plot survival curves by trimming treatment -----------------------------------------------------

# Prepare graphics device
tiff(filename = "Figure 1.tif", width = 2700, height = 2000, units = "px", res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2000, 2700)
pushViewport(viewport(layout = gly))

# Plot unwarmed CN
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CN_NW, bottom = FALSE, left = TRUE, atext = "CN Unwarmed"))
popViewport()

# Plot warmed CN
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CN_W, bottom = TRUE, left = TRUE, atext = "CN Warmed"))
popViewport()

# Plot unwarmed CA
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CA_NW, bottom = FALSE, left = FALSE, atext = "CA Unwarmed"))
popViewport()

# Plot warmed CA
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CA_W, bottom = TRUE, left = FALSE, atext = "CA Warmed"))
popViewport()

# Create legend
grid.text(label = c("Control", "Trim to 10 cm", "Trim to 5 cm", "Trim to Ground"),
          x = rep(0.908, 4), y = seq(0.913, 0.840, length.out = 4),
          hjust = 1, gp = gpar(cex = 0.35))
grid.segments(x0 = rep(0.920, 4), y0 = seq(0.913, 0.840, length.out = 4),
              x1 = rep(0.956, 4), y1 = seq(0.913, 0.840, length.out = 4),
              gp = gpar(lty = c(1, 5, 2, 3), lty = c("solid", "3333", "2323", "1212"), 
                        col = c("black", "purple", "green", "orange")))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

