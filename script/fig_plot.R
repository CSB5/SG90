
figfolder<- figsfoldername
pptfolder <- editablepptfname

template_to_save <- function(ggobject, fname, width, height)
{
  
  #cowplot::save_plot(here(figfolder, fname), ggobject, base_width=width, base_height=height, dpi=600)
  ggsave(here(figsfoldername,fname), plot=ggobject, width=width, height=height, units=c("cm"),limitsize = FALSE)
  
}

template_to_save_sepggplot <- function(plot_names, width, height,site){
  dir.create(here(figfolder, site), showWarnings = FALSE)
  for (i in 1:length(plot_names)) {
    
    cowplot::save_plot(here(figfolder,site, paste0(i, ".png")), plot_names[[i]],
                       base_width=width, base_height=height, dpi=600)
  }
}
makeedtiableppt <- function(imgobj, fname,width, height)
{
  p_dml <- rvg::dml(ggobj = imgobj)
  # Write the document to a file
  officer::read_pptx() %>%
    # add slide ----
  officer::add_slide() %>%
    # specify object and location of object ----
  officer::ph_with(p_dml, officer::ph_location(width = width, height = height)) %>%
    # export slide -----
  base::print(target =
                here(pptfolder, fname))
}

create_dml <- function(plot){
  rvg::dml(ggobj = plot)
}

create_pptx <- function(plot, path, width = width, height = height){
  if (!file.exists(path)) {
    out <- officer::read_pptx()
  }
  else {
    out <- officer::read_pptx(path)
  }
  
  out %>%
    officer::add_slide() %>%
    officer::ph_with(plot, location = officer::ph_location(width = width, height = height)) %>%
    base::print(target = path)
}

makeedtiableppt_manyobj <- function(imgobjlist, fname,width, height){
  pics_dml <- purrr::map(imgobjlist, create_dml)
  
  purrr::map(pics_dml, create_pptx,width = width, height = height, path = here(pptfolder, fname))}
