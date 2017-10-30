load 'load.rb'

model = -1 # -1 means the most recent profile
log_dir = '../LOGS'

pp = ProfilePlots.new(
    'xaxis_by' => 'logP',
    'model' => model,   
    'log_dir' => log_dir,
    'log_mass_frac_ymin' => -6, 
    'log_mass_frac_ymax' => nil,
    'xmin' => 5.8, 'xmax' => 12.5)
    #'xmin' => 0, 'xmax' => 1.4)
    
#load 'extras.rb'
#pp.add_extras
