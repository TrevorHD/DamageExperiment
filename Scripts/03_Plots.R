##### [F1] Plot survival curves by trimming treatment -----------------------------------------------------

# Prepare graphics device
tiff(filename = "Figure 1 NEW.tif", width = 2700, height = 2000, units = "px", res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2000, 2700)
pushViewport(viewport(layout = gly))

# Plot unwarmed CN (year 1)
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 50:1150))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 0, species = "CN", Data_Alt_Y1),
              bottom = FALSE, colm = "left", atext = "CN NW"))
popViewport()

# Plot warmed CN (year 1)
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 50:1150))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 1, species = "CN", Data_Alt_Y1),
              bottom = TRUE, colm = "left", atext = "CN W"))
popViewport()

# Plot unwarmed CA (year 1)
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 1175:2125))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 0, species = "CA", Data_Alt_Y1),
              bottom = FALSE, colm = "middle", atext = "CA NW"))
popViewport()

# Plot warmed CA (year 1)
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 1175:2125))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 1, species = "CA", Data_Alt_Y1),
              bottom = TRUE, colm = "middle", atext = "CA W"))
popViewport()

# Plot unwarmed CA (year 2)
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 2150:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 0, species = "CA", Data_Alt_Y2),
              bottom = FALSE, colm = "right", atext = NULL))
popViewport()

# Plot warmed CA (year 2)
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 2150:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km.curve(comp = "trim", subs = 1, species = "CA", Data_Alt_Y2),
              bottom = TRUE, colm = "right", atext = NULL))
popViewport()

# Create legend
grid.text(label = c("Control", "10 cm", "5 cm", "0 cm"),
          x = rep(0.357, 4), y = seq(0.913, 0.840, length.out = 4),
          hjust = 1, gp = gpar(cex = 0.35))
grid.segments(x0 = rep(0.369, 4), y0 = seq(0.912, 0.837, length.out = 4),
              x1 = rep(0.405, 4), y1 = seq(0.912, 0.837, length.out = 4),
              gp = gpar(lty = c(1, 5, 2, 3), lty = c("solid", "3333", "2323", "1212"), 
                        col = c("green4", "olivedrab3", "gold", "orange")))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()





##### [F2] Plot survival curves by trimming treatment -----------------------------------------------------

# Prepare graphics device
tiff(filename = "Figure 2 NEW.tif", width = 2700, height = 3850, units = "px", res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(3850, 2700)
pushViewport(viewport(layout = gly))

# Plot CN Treatment 1 (year 1)
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 50:1150))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 1, species = "CN", Data_Alt_Y1),
               row = 1, colm = "left", atext = "CN Control"))
popViewport()

# Plot CN Treatment 2 (year 1)
pushViewport(vp = viewport(layout.pos.row = 950:1850, layout.pos.col = 50:1150))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 2, species = "CN", Data_Alt_Y1),
               row = 2, colm = "left", atext = "CN 10 cm"))
popViewport()

# Plot CN Treatment 3 (year 1)
pushViewport(vp = viewport(layout.pos.row = 1875:2775, layout.pos.col = 50:1150))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 3, species = "CN", Data_Alt_Y1),
               row = 3, colm = "left", atext = "CN 5 cm"))
popViewport()

# Plot CN Treatment 4 (year 1)
pushViewport(vp = viewport(layout.pos.row = 2825:3825, layout.pos.col = 50:1150))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 4, species = "CN", Data_Alt_Y1),
               row = 4, colm = "left", atext = "CN 0 cm"))
popViewport()

# Plot CA Treatment 1 (year 1)
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 1175:2125))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 1, species = "CA", Data_Alt_Y1),
               row = 1, colm = "middle", atext = "CA Control"))
popViewport()

# Plot CA Treatment 2 (year 1)
pushViewport(vp = viewport(layout.pos.row = 950:1850, layout.pos.col = 1175:2125))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 2, species = "CA", Data_Alt_Y1),
               row = 2, colm = "middle", atext = "CA 10 cm"))
popViewport()

# Plot CA Treatment 3 (year 1)
pushViewport(vp = viewport(layout.pos.row = 1875:2775, layout.pos.col = 1175:2125))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 3, species = "CA", Data_Alt_Y1),
               row = 3, colm = "middle", atext = "CA 5 cm"))
popViewport()

# Plot CA Treatment 4 (year 1)
pushViewport(vp = viewport(layout.pos.row = 2825:3825, layout.pos.col = 1175:2125))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 4, species = "CA", Data_Alt_Y1),
               row = 4, colm = "middle", atext = "CA 0 cm"))
popViewport()

# Plot CA Treatment 1 (year 2)
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 2150:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 1, species = "CA", Data_Alt_Y2),
               row = 1, colm = "right", atext = NULL))
popViewport()

# Plot CA Treatment 2 (year 2)
pushViewport(vp = viewport(layout.pos.row = 950:1850, layout.pos.col = 2150:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 2, species = "CA", Data_Alt_Y2),
               row = 2, colm = "right", atext = NULL))
popViewport()

# Plot CA Treatment 3 (year 2)
pushViewport(vp = viewport(layout.pos.row = 1875:2775, layout.pos.col = 2150:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 3, species = "CA", Data_Alt_Y2),
               row = 3, colm = "right", atext = NULL))
popViewport()

# Plot CA Treatment 4 (year 2)
pushViewport(vp = viewport(layout.pos.row = 2825:3825, layout.pos.col = 2150:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot2(km.curve(comp = "warm", subs = 4, species = "CA", Data_Alt_Y2),
               row = 4, colm = "right", atext = NULL, blank = TRUE))
popViewport()

# Create legend
grid.text(label = c("Unwarmed", "Warmed"),
          x = rep(0.357, 2), y = c(0.950, 0.937),
          hjust = 1, gp = gpar(cex = 0.35))
grid.segments(x0 = rep(0.369, 2), y0 = c(0.950, 0.937),
              x1 = rep(0.405, 2), y1 = c(0.950, 0.937),
              gp = gpar(lty = c(1, 5, 2, 3), lty = c("solid", "3333"), 
                        col = c("blue", "red")))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

