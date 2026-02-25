


## Colour gradient
grad_cols <- c(
  '#332A82','#3838A2','#306ED9','#55B1AE','#B7BD71','#EDBA5F','#FECC56','#FEE053'
)


white_axis_theme<-  theme(
  # move ticks inside
  axis.ticks.length = unit(0, "cm"),
  
  # move tick labels inside
  axis.text.x = element_text(
    colour = "white",
    margin = margin(t = -14)
  ),
  axis.text.y = element_text(
    colour = "white",
    margin = margin(r = -25)
  ),
  
  # make tick marks white
  axis.ticks = element_line(colour = "white"),
  
  # optional: axis line white too
  axis.line = element_line(colour = "white")
)
