# diff figure 2


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
library(spdep)
library(png)
library(DiagrammeR)
rm(list = setdiff(ls(), lsf.str()))  


# SR bubble ---------------------------------------------------------------

sr <- grViz("
digraph SEM {

  # Node statements
//  subgraph clustercenter{
    node [shape = circle, fontname = 'Helvetica', fontsize=10, style=filled, 
    color=white, peripheries=2, fillcolor=grey90, margin=0]
      SR[label='Species richness \n R&#xb2;=0.65']
//  }
//  subgraph clusterclimate {
    node [shape = rectangle, style=filled, fillcolor = '#71A8F4', 
    color=white, peripheries=2, margin=0.05]
    mat
    pre
//  }
//  subgraph clusterbio {
    node [style=filled, fillcolor = '#9DF471', 
    color=white, peripheries=2, margin=0.01]
    sub_trop_mbf[label='(sub)trop bmf']
    mont_gs[label='mont g&s']
//  }
//  subgraph clusterseason {
    node [style=filled, fillcolor = '#F4F471', 
    color=white, peripheries=2, margin=0.05]
      prs; tra
//  }
//  subgraph clusterEH {
    node [style=filled, fillcolor = '#F471B2', 
    color=white, peripheries=2, margin=0.05]
    soil; area; tri
//  }
    node [style=filled, fillcolor = '#F48771', 
    color=white, peripheries=2]
    mat_ano[label='mat \n anomaly', margin=0.01]

  # Edge statements
  
  edge [color = red, penwidth=1, len=1]
    prs -> SR

  edge [color = red, penwidth=1, style=dashed, len=1]
    tra -> sub_trop_mbf
  
  edge [color = red, penwidth=3, style=normal, len=1, arrowsize=0.7]
    
  edge [color = black, penwidth=1, style=normal, len=1, arrowsize=1]
    {tra tri area} -> SR
    mont_gs -> soil
    {area mat tri}  -> sub_trop_mbf
  
  edge [color = black, penwidth=1, style=dashed, minlen=1]
    {mont_gs} -> SR

  
  edge [color = black, penwidth=3, style=normal, len=1, arrowsize=0.7]
    {mat sub_trop_mbf} -> SR

    area -> soil
  
  edge [color = black, penwidth=7, style=normal, arrowsize=0.4]
    soil -> SR
    pre -> sub_trop_mbf
    
  edge [style=invis]
    pre -> SR
    mat_ano -> SR
    
  graph[layout=twopi, ranksep=2, root=SR, splines=true]
    
}
")


# MRD bubble --------------------------------------------------------------

mrd <- grViz("
digraph SEM2 {

  # Node statements
//  subgraph clustercenter{
    node [shape = circle, fontname = Helvetica, fontsize=10, style=filled, 
    color=white, peripheries=3, fillcolor=grey90, margin=0]
      MRD [label='Diversification \n R&#xb2;=0.56', margin=0]
//  }
//  subgraph clusterclimate {
    node [shape = rectangle, style=filled, fillcolor = '#71A8F4', 
    color=white, peripheries=2, margin=0.05]
    mat
    pre
//  }
//  subgraph clusterbio {
    node [style=filled, fillcolor = '#9DF471', 
    color=white, peripheries=2, margin=0.01]
    sub_trop_mbf[label='(sub)trop bmf']
    mont_gs[label='mont g&s']
//  }
//  subgraph clusterseason {
    node [style=filled, fillcolor = '#F4F471', 
    color=white, peripheries=2, margin=0.05]
      prs; tra
//  }
//  subgraph clusterEH {
    node [style=filled, fillcolor = '#F471B2', 
    color=white, peripheries=2, margin=0.05]
    soil; area; tri
//  }
    node [style=filled, fillcolor = '#F48771', 
    color=white, peripheries=2]
    mat_ano[label='mat \n anomaly', margin=0.01]

  # Edge statements
  
  edge [color = red, penwidth=1, len=1]
    area -> MRD
  
  edge [color = red, penwidth=1, style=dashed, len=1]
    prs -> MRD
    tra -> sub_trop_mbf
  
  edge [color = red, penwidth=3, style=normal, len=1, arrowsize=0.7]
    {pre sub_trop_mbf} -> MRD
    
  edge [color = black, penwidth=1, style=normal, len=1, arrowsize=1]
    mat_ano -> MRD
    mont_gs -> soil
    {area mat tri}  -> sub_trop_mbf
  
  edge [color = black, penwidth=1, style=dashed, minlen=1]
    tri -> MRD
  
  edge [color = black, penwidth=3, style=normal, len=1, arrowsize=0.7]
    soil -> MRD
    area -> soil
  
  edge [color = black, penwidth=7, style=normal, arrowsize=0.4]

    pre -> sub_trop_mbf
    
  edge [style=invis]
  {tra mat mont_gs} -> MRD
    
  graph[layout=twopi, ranksep=2, root=MRD, splines=true]
    
}
")

# 2. Convert to SVG, then save as png
sr = DiagrammeRsvg::export_svg(sr)
sr = charToRaw(sr) # flatten
rsvg::rsvg_svg(sr, "../figures+tables/sr_bubble.svg")
mrd = DiagrammeRsvg::export_svg(mrd)
mrd = charToRaw(mrd) # flatten
rsvg::rsvg_svg(mrd, "../figures+tables/mrd_bubble.svg")
