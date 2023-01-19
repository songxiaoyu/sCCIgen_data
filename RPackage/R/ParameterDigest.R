# Take input

input="../ParameterFile/parameter_file_for_xiaoyu.tsv"

ParameterDigest=function(input){

  para=data.frame(fread(input))
  expression_data_file=para[para$parameters=="expression_data_file", "value"]
  spatial_data_file=para[para$parameters=="spatial_data_file", "value"]

  # active spatial function index
  if (spatial_data_file=="") {
    spatial_sim_idx=1
  } # `PointSimulator_NoData`
  if (spatial_data_file!="" & XX =="XX") {
    spatial_sim_idx=2
  } # `PointSimulator_STData`
  if (spatial_data_file!="" & XX !="XX") {
    spatial_sim_idx=3
  } # Pass existing data


  # pass in parameters for `PointSimulator_NoData`
  if (spatial_sim_idx==1) {
    simulation_seed=as.numeric(para[para$parameters==
                                      "simulation_seed", "value"])
    num_simulated_cells=as.numeric(para[para$parameters==
                                          "num_simulated_cells", "value"])
    num_regions=para[para$parameters==
                       "num_regions", "value"]
    cell_overlap_cutoff=para[para$parameters==
                               "cell_overlap_cutoff", "value"]
    cell_even_distribution=para[para$parameters==
                                  "cell_even_distribution", "value"]

    custom_cell_type_proportions=para[para$parameters==
                                        "custom_cell_type_proportions", "value"]

    if (custom_cell_type_proportions=="n") {
#     XXXXX
    } else {
      cell_type_proportion=para[grep("cell_type_proportion_",
                                     para$parameters), "value"]
    }

    custom_cell_location_interactions=para[para$parameters==
                                             "custom_cell_location_interactions", "value"]
    if (custom_cell_type_proportions=="n") {
      #     XXXXX
    } else {
      #     XXXXX
    }

  } # end of passing  parameters into `PointSimulator_NoData`




}

