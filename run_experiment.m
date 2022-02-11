clear all
close all

format long

% disp('What do you want to do?')
% disp('1 - Solve forward problem')
% disp('2 - Perform simple inversion')
% disp('3 - Compare order of accuracy of forward and adjoint methods')
% %disp('4 - Compute sensitivity matrix')
% disp('4 - Run derivative and adjoint tests')
% choice = input('Enter choice: ');
% while 1
%     if choice > 0 && choice < 5
%         break
%     else
%         choice = input('Try again: ');
%     end
% end
% 
% if choice > 0 && choice < 3
%     disp(' ')
%     disp('Working on a periodic domain [0,40pi] x [0,40pi]. Choose grid size:')
%     disp('1 -  32 x  32')
%     disp('2 -  64 x  64')
%     disp('3 - 128 x 128')
%     choice_grid_size = input('Enter choice: ');
%     while 1
%         if choice_grid_size > 0 && choice_grid_size < 4
%             break
%         else
%             choice_grid_size = input('Try again: ');
%         end
%     end
% else
%     choice_grid_size = 1;
% end

choice = 1;
choice_grid_size = 2;

% =========================================================================

paths('set');

N   = 2^(4+choice_grid_size);
sim = Simulation(N);

result = sim.run(choice);

paths('unset');
