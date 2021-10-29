function new_matrix = searchlight_RSA(simulation, perception, control, cube_size);

new_matrix = zeros(size(simulation));
contrast_dimensions = size(simulation);

for h = 1:contrast_dimensions(3)-cube_size
    for w = 1:contrast_dimensions(2)-cube_size
        for l = 1:contrast_dimensions(1)-cube_size
            simulation_cube = simulation(l:(l+cube_size), w:(w+cube_size), h:(h+cube_size));
            perception_cube = perception(l:(l+cube_size), w:(w+cube_size), h:(h+cube_size));
            control_cube = control(l:(l+cube_size), w:(w+cube_size), h:(h+cube_size));
            
            simulation_list = simulation_cube(:);
            perception_list = perception_cube(:);
            control_list = control_cube(:);
            
            s_p_correlation = corrcoef(simulation_list, perception_list);
            s_p_correlation = s_p_correlation(1, 2);
            
            s_c_correlation = corrcoef(simulation_list, control_list);
            s_c_correlation = s_c_correlation(1, 2);
            
            test = mean(simulation_list);
            
            difference = s_p_correlation - s_c_correlation;
            cube_center_index = [(l+(cube_size/2)) (w+(cube_size/2)) (h+(cube_size/2))];
            cube_center_index = cube_center_index - 0.5;
            new_matrix(cube_center_index(1), cube_center_index(2), cube_center_index(3)) = difference;
            %new_matrix(cube_center_index(1), cube_center_index(2), cube_center_index(3)) = test;
        end
    end
end