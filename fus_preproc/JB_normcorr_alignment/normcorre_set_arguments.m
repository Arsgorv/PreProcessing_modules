function globalArguments = normcorre_set_arguments(globalArguments)
  %{
    Defines various arguments used for various kinds of batch normcorre scripts
    I run. Helps to keep them consistent. 
  %}
if ~isfield(globalArguments,'normcorreFilterType')
  normcorreFilterType = 'defaultSVD';
end
switch globalArguments.normcorreFilterType
  case 'defaultSVD'
    globalArguments.normcorreArguments.max_shift = [3 1]; % don't allow much rigid shifting
    globalArguments.normcorreArguments.max_dev = [3 1];
    globalArguments.normcorreArguments.mot_uf = 4;
    globalArguments.normcorreArguments.vary_grid_size = 1;
    if ~globalArguments.normcorreArguments.vary_grid_size
      error('Don''t it''s a bad idea'); % I have this here because if I change my mind I want to avoid forgetting I reset grid_size each loop. But I can't imagine changing my mind right now. 
      %grid_size = [size(b0_in,1),6]; %%% commented out because I prefer to set this within each loop
    else

    end
    globalArguments.normcorreArguments.overlap_pre = [48 24]; %%% So I intentionally overlap my patches a large amount because I want to keep the corrections mostly rigid. 
    globalArguments.normcorreArguments.overlap_post = [48 24]; %%% pre and post have to do with pre and post upsampling. Despite that, they are in the same units, physics-wise. 
    globalArguments.normcorreArguments.shifts_method = 'linear';
    globalArguments.normcorreArguments.init_batch = 2; %%% just use two frames. Image is clean, and this way I don't take mean over moving images. Won't let me use one frame, but 2 should be fine. 
    globalArguments.normcorreArguments.min_diff = [16,16];
    globalArguments.normcorreArguments.correct_bidir = false;
  case {'tissue','keep_all'}
    globalArguments.normcorreArguments.max_shift = [3 1]; % don't allow much rigid shifting
    globalArguments.normcorreArguments.max_dev = [3 1];
    globalArguments.normcorreArguments.mot_uf = 4;
    globalArguments.normcorreArguments.vary_grid_size = 1;
    if ~globalArguments.normcorreArguments.vary_grid_size
      %%% vary grid size for each slice, somehow. Likely, make it relative to the size. 
    else
      globalArguments.normcorreArguments.grid_size = [4,6]; %%% still might prefer to set this in a size-dependent manner because pixels correspond to a set physical distance... If I do that I should set "vary_grid_size" to 1 and implement it in-loop. 
    end
    globalArguments.normcorreArguments.overlap_pre = [12 12]; %%% Testing less overlap for the tissue NoRMcorre
    globalArguments.normcorreArguments.overlap_post = [12 12]; %%% pre and post have to do with pre and post upsampling. Despite that, they are in the same units, physics-wise. 
    globalArguments.normcorreArguments.shifts_method = 'linear';
    globalArguments.normcorreArguments.init_batch = 2; %%% just use two frames. Image is clean, and this way I don't take mean over moving images. Won't let me use one frame, but 2 should be fine. 
    globalArguments.normcorreArguments.min_diff = [16,16];
    globalArguments.normcorreArguments.correct_bidir = false;
end
	%% both the following are irrelevant to the new ferrets. 
	  %% set a constant
	    globalArguments.normcorreArguments.ferretDependent.minFramesForGoodRecording = 9; % less than 9, recording is unusable. %%% this is irrelevant to the new ferrets.
	  %% set a very important argument!!!
	    globalArguments.normcorreArguments.ferretDependent.framesIncluded = 2:9; % decides which frames will be included in everything going forward. At time of writing, 2:9 is the best option. 
end
  
  
  
  
  
  
  
