%binary_file_sandbox.m
%
% Purpose: read in Sounder binary files and understand the format
%
% Use data from binary file to created frequency domain (FD) file.  Use FD
% file to create time domain (TD) file.
%
% ul_data_b_16QAM_52_64_10_1_1_AB_0.bin % Binary data
% ul_data_f_16QAM_52_64_10_1_1_AB_0.bin % FD constellations
% ul_data_t_16QAM_52_64_10_1_1_AB_0.bin % TD over-the-air
%
% Notes from https://wiki.renew-wireless.org/en/architectures/softwarearchitectureoverview:
%  16QAM
%  52 active subcarriers in OFDM symbol
%  64 total subcarriers in OFDM symbol
%  10 OFDM symbols per slot
%  1 uplink data slot
%  1 frame
%  A and B channels
%  0 is the client index
%
% Notes from rl_ofdm_mimo.m:
% SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
% SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
% SC_IND_DATA_PILOT       = [2:27 39:64]'; % Nothing in the center and 5 or 6 empty bins on the sides of the band
%

%% Binary Data
% 8320 bits were read
filename = 'ul_data_b_16QAM_52_64_10_1_1_AB_0.bin';
file_id = fopen( filename );
temp_data = fread( file_id, 'ubit1' );
fclose( file_id );

% Remove even nibbles
stream1 = temp_data( 1 : 8 : end );
stream2 = temp_data( 2 : 8 : end );
stream3 = temp_data( 3 : 8 : end );
stream4 = temp_data( 4 : 8 : end );
stream5 = temp_data( 5 : 8 : end );
stream6 = temp_data( 6 : 8 : end );
stream7 = temp_data( 7 : 8 : end );
stream8 = temp_data( 8 : 8 : end );
intermediate_data = upsample( stream1, 4 ) + ...
    circshift( upsample( stream2, 4 ), 1 ) + ...
    circshift( upsample( stream3, 4 ), 2 ) + ...
    circshift( upsample( stream4, 4 ), 3 );
zero_data = upsample( stream5, 4 ) + ...
    circshift( upsample( stream6, 4 ), 1 ) + ...
    circshift( upsample( stream7, 4 ), 2 ) + ...
    circshift( upsample( stream8, 4 ), 3 );

% The first 48 symbols are valid and the last 4 are invalid
binary_data = nan( 4 * 48 * 10 * 2, 1 );
payload_length = 4 * 48;
frame_length = 4 * 52;
for idx = 1 : 20
    binary_data( ( idx - 1 ) * payload_length + ( 1 : payload_length ) ) = ...
        intermediate_data( ( idx - 1 ) * frame_length + ( 1 : payload_length ) );
end
disp( [ 'Number of bits: ' num2str( numel( binary_data ) ) ] )
% NOTE: The data comes in frames of 192 bits

% % % Plot binary data
% % figure( 101 ), clf
% % plot( binary_data )
% % axis( [ 1 numel( binary_data ) -1 2 ] )
% % grid on

% Convert data to 16QAM
M = 16; % modulation order
map_vect = [ 4; 12; 8; 0; 6; 14; 10; 2; 7; 15; 11; 3; 5; 13; 9; 1 ]; % This is a different version of gray coding than MATLAB's default
x = qammod( binary_data, M, map_vect, 'InputType', 'bit', 'UnitAveragePower',true ); % transmit data
% % z = qamdemod( x, M, 'gray', 'OutputType', 'bit', 'UnitAveragePower',true); % binary data
figure( 102 ), clf
plot( real( x ), imag( x ), 'o' )
axis equal
grid on

% Plot modulated data (i.e., FD)
figure( 103 ), clf
symbol_data_length = 48;
subplot( 211 ), plot( real( x( 1 : symbol_data_length ) ), 'o-' )
grid on
subplot( 212 ), plot( imag( x( 1 : symbol_data_length ) ), 'o-' )
grid on


%% FD
% This generates 20 frames of data (64 symbols per frame).
% There are 4 pilot symbols per frame and 48 data symbols (i.e. 192 bits)
filename = 'ul_data_f_16QAM_52_64_10_1_1_AB_0.bin';
file_id = fopen( filename );
freq_data = fread( file_id, 'single' );
fclose( file_id );
i_freq_data = freq_data( 1 : 2 : end );
q_freq_data = freq_data( 2 : 2 : end );
fd_data = i_freq_data + 1j * q_freq_data; % complex data
disp( [ 'Number of symbols: ' num2str( numel( fd_data ) ) ] )

% Duplicate FD data from binary data
%  Map symbols to subcarriers and insert pilots
data_symbols = [ 7 : 11, 13 : 25, 27 : 32, 34 : 39, 41 : 53, 55 : 59 ];
pilot_symbols = [ 12, 26, 40, 54 ];
pilots = [ 1, 1, -1, 1 ];
fd_data_check = zeros( 1280, 1 );
for idx = 1 : 20 % for each frame
    fd_data_check( ( idx - 1 ) * 64 + data_symbols ) = x( ( idx - 1 ) * 48 + ( 1 : 48 ) );
    fd_data_check( ( idx - 1 ) * 64 + pilot_symbols ) = pilots;
end
fd_data_error = sum( abs( fd_data_check - fd_data ) ); % This value is very small

% Plot constellation
figure( 201 )
plot( i_freq_data, q_freq_data, 'o' )
axis equal
grid on

% Plot FD Data
figure( 202 )
subplot( 211 ), plot( i_freq_data, 'o-' )
grid on
subplot( 212 ), plot( q_freq_data, 'o-' )
grid on

% Plot first OFDM symbol (data only)
figure( 203 )
symbol_length = 64;
data_symbols = [ 7 : 11, 13 : 25, 27 : 32, 34 : 39, 41 : 53, 55 : 59 ];
subplot( 211 ), plot( i_freq_data( data_symbols ), 'o-' )
grid on
subplot( 212 ), plot( q_freq_data( data_symbols ), 'o-' )
grid on
mod_data = [];
for idx = 1 : 20
    mod_data = [ mod_data; i_freq_data( ( idx - 1 ) * 64 + data_symbols ) + ...
        1j * q_freq_data( ( idx - 1 ) * 64 + data_symbols ) ];
end
z = qamdemod( mod_data, M, map_vect, 'OutputType', 'bit', 'UnitAveragePower', true); % binary data
fd_data_error2 = sum( abs( binary_data - z ) ); % This checks out as expected

% Convert the first symbol to the TD
% This is correct except for the scaling
first_symbol_FD = fftshift( fd_data( 1 : symbol_length ) );
first_symbol_TD = ifft( first_symbol_FD ) * 2^15; % * sqrt( symbol_length );
figure( 204 )
subplot( 211 ), plot( real( first_symbol_TD ), 'o-' )
grid on
subplot( 212 ), plot( imag( first_symbol_TD ), 'o-' )
grid on

%% TD
filename = 'ul_data_t_16QAM_52_64_10_1_1_AB_0.bin';
file_id = fopen( filename );
time_data = fread( file_id, 'int16' );
fclose( file_id );
i_time_data = time_data( 1 : 2 : end );
q_time_data = time_data( 2 : 2 : end );
td_data = i_time_data + 1j * q_time_data;

symbol_length = 64;
cp_length = 16;
num_symbols = 10;
pad_length = 160;
num_clients = 2;

% Duplicate TD data from FD data
%  NOTE: 160 samples of pad flank the 10 OFDM symbols for each client
td_data_check = zeros( num_clients * ( 2 * pad_length + num_symbols * ( symbol_length + cp_length ) ), 1 );
vect_idx = 1;
sym_idx = 1;
for client_idx = 1 : num_clients
    vect_idx = vect_idx + pad_length; % advance pointer for pad
    for idx = 1 : num_symbols
        curr_fd_sym = fd_data( sym_idx + ( 0 : symbol_length - 1 ) );
        sym_idx = sym_idx + symbol_length;
        curr_td_sym = ifft( fftshift( curr_fd_sym ) ) * 2^15; % The shift moves the zero-frequency bin to the beginning of the vector
        % Write the CP
        td_data_check( vect_idx + ( 0 : cp_length - 1 ) ) = curr_td_sym( end - cp_length + 1 : end );
        vect_idx = vect_idx + cp_length;
        % Write the symbol
        td_data_check( vect_idx + ( 0 : symbol_length - 1 ) ) = curr_td_sym; 
        vect_idx = vect_idx + symbol_length;
    end
    vect_idx = vect_idx + pad_length; % advance pointer for pad
end
figure( 300 ), clf
subplot( 211 ), plot( real( td_data ) )
hold on, grid on
plot( real( td_data_check ) )
hold off
subplot( 212 ), plot( imag( td_data ) )
hold on, grid on
plot( imag( td_data_check ) )
hold off



offset1 = pad_length;
offset2 = 3 * pad_length + num_symbols * ( symbol_length + cp_length );

% Plot TD symbols
figure( 301 ), clf
subplot( 211 ), plot( i_time_data, 'o-' )
grid on
hold on % plot CP in red
for idx = 1 : num_symbols
    plot_indices = offset1 + ( idx * symbol_length ) + ( idx - 1 ) * cp_length + ( 1 : cp_length );
    plot( plot_indices, i_time_data( plot_indices ), 'ro-' )
end
hold off
subplot( 212 ), plot( q_time_data, 'o-' )
grid on
hold on % plot CP in red
for idx = 1 : num_symbols
    plot_indices = offset1 + ( idx * symbol_length ) + ( idx - 1 ) * cp_length +  ( 1 : cp_length );
    plot( plot_indices, q_time_data( plot_indices ), 'ro-' )
end
hold off

%% Determine gain value between FD conversion and TD data
symbol_length = 64;
cp_length = 16;
offset1 = 160;
first_symbol_FD = fftshift( i_freq_data( 1 : symbol_length ) + 1j * q_freq_data( 1 : symbol_length ) );
first_symbol_TD = ifft( first_symbol_FD ) * sqrt( symbol_length );
idx = 1;
plot_indices = offset1 + ( idx - 1 ) * ( symbol_length + cp_length ) + cp_length + ( 1 : symbol_length );
TD_check = i_time_data( plot_indices ) + 1j * q_time_data( plot_indices );
real_scale_vect = real( TD_check ) ./ real( first_symbol_TD );
imag_scale_vect = imag( TD_check ) ./ imag( first_symbol_TD );
scale_factor = mean( [ real_scale_vect; imag_scale_vect ] );
TD_test = round( real( first_symbol_TD ) * scale_factor ) + ...
    1j * round( imag( first_symbol_TD ) * scale_factor );
error_vect = TD_test - TD_check;






