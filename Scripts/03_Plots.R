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
print(km.plot(km.curve(comp = "trim", subs = 0, species = "CN", Data_Alt_Y1),
              bottom = FALSE, left = TRUE, atext = "CN Unwarmed"))
popViewport()

# Plot warmed CN
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 1, species = "CN", Data_Alt_Y1),
              bottom = TRUE, left = TRUE, atext = "CN Warmed"))
popViewport()

# Plot unwarmed CA
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 0, species = "CA", Data_Alt),
              bottom = FALSE, left = FALSE, atext = "CA Unwarmed"))
popViewport()

# Plot warmed CA
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 1, species = "CA", Data_Alt),
              bottom = TRUE, left = FALSE, atext = "CA Warmed"))
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





##### [F2] Plot survival curves by trimming treatment -----------------------------------------------------

# Prepare graphics device
tiff(filename = "Figure 2.tif", width = 2700, height = 3900, units = "px", res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(4000, 2700)
pushViewport(viewport(layout = gly))

# Plot CN Treatment 1
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 1, species = "CN", Data_Alt_Y1),
               row = 1, left = TRUE, atext = "CN Control"))
popViewport()

# Plot CN Treatment 2
pushViewport(vp = viewport(layout.pos.row = 960:1860, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 2, species = "CN", Data_Alt_Y1),
               row = 2, left = TRUE, atext = "CN Trim to 10 cm"))
popViewport()

# Plot CN Treatment 3
pushViewport(vp = viewport(layout.pos.row = 1910:2810, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 3, species = "CN", Data_Alt_Y1),
               row = 3, left = TRUE, atext = "CN Trim to 5 cm"))
popViewport()

# Plot CN Treatment 4
pushViewport(vp = viewport(layout.pos.row = 2875:3875, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 4, species = "CN", Data_Alt_Y1),
               row = 4, left = TRUE, atext = "CN Trim to Ground"))
popViewport()

# Plot CA Treatment 1
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 1, species = "CA", Data_Alt),
               row = 1, left = FALSE, atext = "CA Control"))
popViewport()

# Plot CA Treatment 2
pushViewport(vp = viewport(layout.pos.row = 960:1860, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 2, species = "CA", Data_Alt),
               row = 2, left = FALSE, atext = "CA Trim to 10 cm"))
popViewport()

# Plot CA Treatment 3
pushViewport(vp = viewport(layout.pos.row = 1910:2810, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 3, species = "CA", Data_Alt),
               row = 3, left = FALSE, atext = "CA Trim to 5 cm"))
popViewport()

# Plot CA Treatment 4
pushViewport(vp = viewport(layout.pos.row = 2875:3875, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 4, species = "CA", Data_Alt),
               row = 4, left = FALSE, atext = "CA Trim to Ground"))
popViewport()

# Create legend
grid.text(label = c("Unwarmed", "Warmed"),
          x = rep(0.908, 2), y = c(0.950, 0.937),
          hjust = 1, gp = gpar(cex = 0.35))
grid.segments(x0 = rep(0.920, 2), y0 = c(0.950, 0.937),
              x1 = rep(0.956, 2), y1 = c(0.950, 0.937),
              gp = gpar(lty = c(1, 5, 2, 3), lty = c("solid", "3333"), 
                        col = c("blue", "red")))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

