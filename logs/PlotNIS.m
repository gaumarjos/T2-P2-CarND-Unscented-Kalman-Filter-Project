function PlotNIS(laser_nis, radar_nis)
  subplot(1,2,1)
  n_laser = length(laser_nis);
  plot(laser_nis);
  hold on;
  plot([0 n_laser], [5.991 5.991]);
  ylim([0 5.991*2])
  
  subplot(1,2,2)
  n_radar = length(radar_nis);
  plot(radar_nis);
  hold on;
  plot([0 n_radar], [7.815 7.815]);
  ylim([0 7.815*2])