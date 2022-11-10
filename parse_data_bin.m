clear all;
close all;

%Read in bin data as uint8 placed in double array
bin_file_name = 'binary_data.bin';
fileID = fopen(['binary_data.bin']);
tx_data = fread(fileID);

num_frames = 2000;
num_ul_slots = 1;
ul_slots = 1;
num_cl_sdr_ch = 2;
sym_per_slot = 10;
fft_size = 64;
cp_size = 16;


ofdm_tx_zero_prefix = 160;
ofdm_tx_zero_postfix = 160;
prefix_zpad = zeros(1, ofdm_tx_zero_prefix_arr);
postfix_zpad = zeros(1, ofdm_tx_zero_postfix_arr);


ue_tx_gain_a = 65;
ue_rx_gain_b = 81;
ue_tx_gain_a = 65;
ue_rx_gain_b = 81;


lts_data_ind = [ 6,  7,  8,  9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,...
                23, 24, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 40, 41,...
                42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58];

lts_pilot_ind = [11, 25, 39, 53];
lts_pilot_val = [1.0, 1.0, -1.0, 1.0];


%% setup modulation scheme/coding 16QAM
scale1 = 1/sqrt(10);
mod_16qam = [-3*scale1 -scale1 scale1 3*scale1];
qam16_table = zeros(2, 4);
for idx = 1:16
	qam16_table(1,idx) = mod_16qam(floor((idx-1)/4)+1);
	qam16_table(2,idx) = mod_16qam(mod((idx-1), 4)+1);
end

% |  03  07  11  15  |
% |  01  06  10  14  |
% |  02  05  09  13  |  
% |  00  04  08  12  |

% |  0011   0111   1011   1111   |  
% |  0001   0110   1010   1110   | 
% |  0010   0101   1001   1101   | 
% |  0000   0100   1000   1100   | 

%figure()
%hold on;
%for ii = 1:16
%
%	x_re = qam16_table(1, ii);
%	y_im = qam16_table(2, ii);
%	plot(x_re, y_im, 'r*');
%end



% I guess the sequence is repeated for each frame?
% Number of bytes in bin file calculation: 
% total_bytes = uplink_slots * client_sdr_channels * symbols_per_slot * data_subcarriers
% total_bytes =    1         *          2          *       10         *      52
chunk52 = -1; % Start at 0 in loop, I hate 1 based indexing!
data_freq_dom = single([]);
data_time_dom = single([]);
for f = 1:num_frames
  for ul = 1:num_ul_slots
  	for ch = 1:num_cl_sdr_ch
  	  for s = 1:sym_per_slot

        chunk52 = chunk52 + 1;
  	  	curr_data = tx_data(52*chunk52+1:52*chunk52+52);

  	  	data_time_dom = [prefix_zpad data_time_dom];
    
        %----------------------------------------------------------
        %  BEGIN
        %  Sounder commlib.cc modulate block QAM16
        %----------------------------------------------------------
        
        % Sounder uses complex float (single precision 4byte, 32 bit, 7 decimal)
        % std::complex<float>(qam16_table[0][in[i]], qam16_table[1][in[i]]);
        
        % C++ float = single Matlab
        % C++ cfloat = single Matlab with complex values
        mod_data = single(zeros(1, 52));
        for ii = 1:52
        
        	% Map uint16 in matlab double to QAM16 symbols type complex float (single)
        	mod_data(ii) = single(qam16_table(1, curr_data(ii)) + ...
        		                i*qam16_table(2, curr_data(ii)));
        
        end
        
        %----------------------------------------------------------
        %  END
        %  Sounder commlib.cc modulate block QAM16
        %----------------------------------------------------------
        
        % Fill data subcarrier indices
        ofdm_sym = single(zeros(1, 64));
        for ii = 1:48
          ofdm_sym(lts_data_ind(ii)) = mod_data(ii);
        end

        % Fill pilot subcarrier indices
        for ii = 1:4
          ofdm_sym(lts_pilot_ind(ii)) = lts_pilot_val(ii);
        end

        % Append FD symbol to array for final output
        data_freq_dom = [data_freq_dom ofdm_sym];

        % std::vector<std::complex<float>> CommsLib::IFFT(
        %   const std::vector<std::complex<float>>& in, int fftSize, float scale,
        %   bool normalize, bool fft_shift) {

        % Sounder ifft takes cfloat and returns cfloat
        % Sounder uses fft_shifting -> fftshift(*)
        % Sounder scales by factor (as float -> single in matlab): 1/fft_size
        % Sounder does NOT normalize
        tx_sym = fftshift(single(1/fft_size)*single(ifft(ofdm_sym, fft_size)));

        % Append the cp length tail of tx syms to beginning of array
        tx_sym = [tx_sym(end-(cp_size-1):end) tx_sym];

        % Append Time Domain Symbol to output array
        data_time_dom = [data_time_dom tx_sym];

        % Convert time domain cfloat -> c

      end % Symbol per slot

      % Insert

    end % sdr channels
  end % ul slots
end % frames

%  static constexpr size_t lts_data_ind[48] = {
%      6,  7,  8,  9,  10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
%      23, 24, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 40, 41,
%      42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58};


%  data_ind_ = CommsLib::getDataSc(fft_size_, symbol_data_subcarrier_num_);

%  std::vector<size_t> data_sc;
%  if (fftSize == Consts::kFftSize_80211) {
%    // We follow 802.11 PHY format here
%    data_sc.assign(Consts::lts_data_ind,
%                   Consts::lts_data_ind + Consts::kNumDataSubcarriers_80211);
%  } else {  // Allocate the center subcarriers as data
%    size_t start_sc = (fftSize - DataScNum) / 2;
%    size_t stop_sc = start_sc + DataScNum;
%    for (size_t i = start_sc; i < stop_sc; i++)
%      if ((i - start_sc) % kPilotSubcarrierSpacing != PilotScOffset)
%        data_sc.push_back(i);
%  }

% Debug/Parse byte file in blocks of 52 bytes
% symii = zeros(1, length(tx_data));
% figure();
% hold on;
% for ii = 1:52:length(tx_data)-51
%   symii = tx_data(ii:ii+51);
%   plot(symii);
% end





