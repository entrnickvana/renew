clear all;
close all;

%Read in bin data as uint8 placed in double array
bin_file_name = 'binary_data.bin';
fileID = fopen(['binary_data.bin']);
tx_data = fread(fileID);

fileID_FD_REF = fopen(['ul_data_f_16QAM_52_64_10_1_1_AB_0.bin']);
data_freq_dom_ref = fread(fileID_FD_REF, 'float');
fclose(fileID_FD_REF);

num_frames = 1;
num_ul_slots = 1;
ul_slots = 1;
num_cl_sdr_ch = 2;
sym_per_slot = 10;
fft_size = 64;
cp_size = 16;


tx_gain_a=81;

ofdm_tx_zero_prefix = 160;
ofdm_tx_zero_postfix = 160;
prefix_zpad = single(zeros(1, ofdm_tx_zero_prefix));
postfix_zpad = single(zeros(1, ofdm_tx_zero_postfix));


ue_tx_gain_a = 65;
ue_rx_gain_b = 81;
ue_tx_gain_a = 65;
ue_rx_gain_b = 81;


% 48 data subcarrier indices
lts_data_ind_cpp =...
                [ 6,  7,  8,  9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,...
                23, 24, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 40, 41,...
                42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58];

lts_data_ind = [ 7,  8,  9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,...
                24, 25, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 41, 42,...
                43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 59];


% Used to (among other thinsg), calculate tx_scale
% C++ from sounder is type cfloat, single is 2 byte complex float
lts_seq = ... 
      single([(0 + i*0)  (0 + i*0)  (0 + i*0)  (0 + i*0)  (0 + i*0)  (0 + i*0)  (1 + i*0)  (1 + i*0)...
              (-1 + i*0) (-1 + i*0) (1 + i*0)  (1 + i*0)  (-1 + i*0) (1 + i*0)  (-1 + i*0) (1 + i*0)...
              (1 + i*0)  (1 + i*0)  (1 + i*0)  (1 + i*0)  (1 + i*0)  (-1 + i*0) (-1 + i*0) (1 + i*0)...
              (1 + i*0)  (-1 + i*0) (1 + i*0)  (-1 + i*0) (1 + i*0)  (1 + i*0)  (1 + i*0)  (1 + i*0)...
              (0 + i*0)  (1 + i*0)  (-1 + i*0) (-1 + i*0) (1 + i*0)  (1 + i*0)  (-1 + i*0) (1 + i*0)...
              (-1 + i*0) (1 + i*0)  (-1 + i*0) (-1 + i*0) (-1 + i*0) (-1 + i*0) (-1 + i*0) (1 + i*0)...
              (1 + i*0)  (-1 + i*0) (-1 + i*0) (1 + i*0)  (-1 + i*0) (1 + i*0)  (-1 + i*0) (1 + i*0)...
              (1 + i*0)  (1 + i*0)  (1 + i*0)  (0 + i*0)  (0 + i*0)  (0 + i*0)  (0 + i*0)  (0 + i*0)]);                


fileID_BIN_OUT = fopen('bin_out.bin', 'w');
fileID_TD_OUT = fopen('td_out.bin', 'w');
fileID_FD_OUT = fopen('fd_out.bin', 'w');

%lts_seq_cint16 = cast(lts_seq*32768, "int16");

% This section reproduces the variable tx_scale
lts_seq_t = (1/fft_size)*fftshift(ifft(lts_seq, fft_size));
lts_seq_cint16 = cast(lts_seq_t*32768, "int16");

max_mag_pilot = max(abs(lts_seq_t));

%lts_seq_cint16_norm_scaled_6dB = lts_seq_cint16/(4*max_mag_pilot)

tx_scale = 1/(4*max_mag_pilot);

%figure()
%stem(lts_seq_cint16_norm_scaled_6dB);
%hold on;
%lts_seq_cint16_norm_scaled_6dB = lts_seq_cint16/(cast(4*max_mag_pilot,"int16"))
%stem(lts_seq_cint16_norm_scaled_6dB);
%hold off;


% 4 Pilot subcarrier indices
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

% 4  6   7  5
% 12 14 15 13
% 8  10 11  9
% 0  2   3  1

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

  	  %% prepend padding to time domain
  	  data_time_dom = [prefix_zpad data_time_dom];

  	  for s = 1:sym_per_slot

        % increment chunks of 52,
        chunk52 = chunk52 + 1;
  	  	curr_data = tx_data(52*chunk52+1:52*chunk52+52);

        fwrite(fileID_BIN_OUT, curr_data, 'uint8');  	  	

        %----------------------------------------------------------
        %  BEGIN
        %  Sounder commlib.cc modulate block QAM16
        %----------------------------------------------------------
        
        % Sounder uses complex float (single precision 4byte, 32 bit, 7 decimal)
        % std::complex<float>(qam16_table[0][in[i]], qam16_table[1][in[i]]);
        
        % C++ float = single Matlab
        % C++ cfloat = single Matlab with complex values
        mod_data = single(zeros(1, 48));
        for ii = 1:48

        	% Map uint16 in matlab double to QAM16 symbols type complex float (single)
        	mod_data(ii) = single(qam16_table(1, curr_data(ii)+1) + i*qam16_table(2, curr_data(ii)+1));
        
        end
        
        %----------------------------------------------------------
        %  END
        %  Sounder commlib.cc modulate block QAM16
        %----------------------------------------------------------
        
        % Fill data subcarrier indices
        ofdm_sym = single(zeros(1, 64));
        for ii = 1:48
          % ofdm_sym(lts_data_ind(ii)+1) = mod_data(ii);        	
          ofdm_sym(lts_data_ind(ii)) = mod_data(ii);
        end

        % Fill pilot subcarrier indices
        for ii = 1:4
          %ofdm_sym(lts_pilot_ind(ii)+1) = lts_pilot_val(ii);        	
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
        % Sounder does NOT normalize in IFFT
        tx_sym = fftshift(single(1/fft_size)*single(ifft(ofdm_sym, fft_size)));

        % Append the cp length tail of tx syms to beginning of array
        tx_sym = [tx_sym(end-(cp_size-1):end) tx_sym];

        % Append Time Domain Symbol to output array
        data_time_dom = [data_time_dom tx_sym];

        % Convert time domain cfloat -> c

      end % Symbol per slot

      % Append padding to time domain
      data_time_dom = [data_time_dom postfix_zpad];
      data_time_dom_scaled_6dB = (1/(4*max_mag_pilot))*data_time_dom;

            % Check for saturation in completed time domain data
      if(single(max(single(data_time_dom_scaled_6dB))) > single(1.0))
        display('Error!!!, time domain data saturated above 1.0');
      end

      % data_time_dom_scaled_6dB_cint16 = cast(data_time_dom_scaled_6dB, "int16")
      data_time_dom_scaled_6dB_cint16 = cast(data_time_dom_scaled_6dB*32768, "int16");

      % normalize tx time domain against max pilot

      %for (size_t id = 0; id < data_time_dom.size(); id++) {
      %  data_time_dom.at(id) *= cfg_->tx_scale();
      %  if (data_time_dom.at(id).real() > 1.f ||
      %      data_time_dom.at(id).imag() > 1.f) {
      %    std::printf(
      %        "Saturation detected in frame %zu slot %zu channel %zu "
      %        "sample %zu\n",
      %        f, u, h, id);
      %  }
      %}

      %% https://www.mathworks.com/help/matlab/ref/int16.html       
      %% Variables in MATLAB® of data type (class) int16 are stored as 2-byte (16-bit)
      %% signed integers. For example:

      %cint16 -> cfloat
      %%https://www.mathworks.com/matlabcentral/answers/363182-converting-complex-int16-matrix-to-32-bit-floating-points-or-double      
      
            
      % std::vector<std::complex<int16_t>> Utils::cfloat_to_cint16(
      %     const std::vector<std::complex<float>>& in) {
      %   size_t len = in.size();
      %   std::vector<std::complex<int16_t>> out(len, 0);
      %   for (size_t i = 0; i < len; i++)
      %     out.at(i) = std::complex<int16_t>((int16_t)(in.at(i).real() * 32768),
      %                                       (int16_t)(in.at(i).imag() * 32768));
      %   return out;
      %}

      %% Write frequency domain data (matlab single 8-byte (4 real, 4 imag))
      fwrite(fileID_FD_OUT, data_freq_dom, 'float'); 
      fwrite(fileID_TD_OUT, data_time_dom_scaled_6dB_cint16, 'uint16');      

    end % sdr channels
  end % ul slots
end % frames

fclose(fileID_BIN_OUT);
fclose(fileID_FD_OUT);
fclose(fileID_TD_OUT);

%fileID_BIN_OUT_COMP = fopen(['bin_out.bin']);
%tx_data_compare = fread(fileID_BIN_OUT_COMP, 'uint8');
fileID_BIN_OUT_COMP = fopen(['ul_data_b_16QAM_52_64_10_1_1_AB_0.bin']);
tx_data_compare = fread(fileID_BIN_OUT_COMP, 'uint8');

fileID_FD_OUT_ORIG = fopen(['fd_out.bin']);
data_freq_dom_orig = fread(fileID_FD_OUT_ORIG, 'float');
fileID_FD_OUT_COMP = fopen(['ul_data_f_16QAM_52_64_10_1_1_AB_0.bin']);
data_freq_dom_compare = fread(fileID_FD_OUT_COMP, 'float');

fileID_TD_OUT_COMP = fopen(['td_out.bin']);
data_time_dom_compare = fread(fileID_TD_OUT_COMP, 'int16');

fclose(fileID_BIN_OUT_COMP);
fclose(fileID_FD_OUT_COMP);
fclose(fileID_TD_OUT_COMP);
fclose(fileID_FD_OUT_ORIG);

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





% Constants from sounder


%class Consts {
% public:
%  static constexpr size_t kFftSize_80211 = 64;
%  static constexpr size_t kNumMappedSubcarriers_80211 = 52;
%  static constexpr size_t kNumDataSubcarriers_80211 = 48;
%  static constexpr size_t kNumPilotSubcarriers_80211 = 4;
%  static constexpr size_t kNumNullSubcarriers_80211 = 12;
%  // Define freq-domain STS according to
%  // https://standards.ieee.org/standard/802_11a-1999.html
%  static constexpr std::complex<float> sts_seq[64] = {
%      {0, 0},   {0, 0}, {0, 0}, {0, 0}, {0, 0},   {0, 0}, {0, 0}, {0, 0},
%      {1, 1},   {0, 0}, {0, 0}, {0, 0}, {-1, -1}, {0, 0}, {0, 0}, {0, 0},
%      {1, 1},   {0, 0}, {0, 0}, {0, 0}, {-1, -1}, {0, 0}, {0, 0}, {0, 0},
%      {-1, -1}, {0, 0}, {0, 0}, {0, 0}, {1, 1},   {0, 0}, {0, 0}, {0, 0},
%      {0, 0},   {0, 0}, {0, 0}, {0, 0}, {-1, -1}, {0, 0}, {0, 0}, {0, 0},
%      {-1, -1}, {0, 0}, {0, 0}, {0, 0}, {1, 1},   {0, 0}, {0, 0}, {0, 0},
%      {1, 1},   {0, 0}, {0, 0}, {0, 0}, {1, 1},   {0, 0}, {0, 0}, {0, 0},
%      {1, 1},   {0, 0}, {0, 0}, {0, 0}, {0, 0},   {0, 0}, {0, 0}, {0, 0}};
%
%  // Define freq-domain LTS according to
%  // https://standards.ieee.org/standard/802_11a-1999.html
%  static constexpr std::complex<float> lts_seq[64] = {
%      {0, 0},  {0, 0},  {0, 0},  {0, 0},  {0, 0},  {0, 0},  {1, 0},  {1, 0},
%      {-1, 0}, {-1, 0}, {1, 0},  {1, 0},  {-1, 0}, {1, 0},  {-1, 0}, {1, 0},
%      {1, 0},  {1, 0},  {1, 0},  {1, 0},  {1, 0},  {-1, 0}, {-1, 0}, {1, 0},
%      {1, 0},  {-1, 0}, {1, 0},  {-1, 0}, {1, 0},  {1, 0},  {1, 0},  {1, 0},
%      {0, 0},  {1, 0},  {-1, 0}, {-1, 0}, {1, 0},  {1, 0},  {-1, 0}, {1, 0},
%      {-1, 0}, {1, 0},  {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {1, 0},
%      {1, 0},  {-1, 0}, {-1, 0}, {1, 0},  {-1, 0}, {1, 0},  {-1, 0}, {1, 0},
%      {1, 0},  {1, 0},  {1, 0},  {0, 0},  {0, 0},  {0, 0},  {0, 0},  {0, 0}};
%
%  static constexpr size_t lts_data_ind[48] = {
%      6,  7,  8,  9,  10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
%      23, 24, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 40, 41,
%      42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58};
%
%  static constexpr size_t lts_pilot_ind[4] = {11, 25, 39, 53};
%  static constexpr std::complex<float> lts_pilot_val[4] = {1.0, 1.0, -1.0, 1.0};
%  static constexpr size_t lts_null_ind[12] = {0,  1,  2,  3,  4,  5,
%                                              32, 59, 60, 61, 62, 63};
%
%  // prime numbers [1,2048]
%  static constexpr size_t prime[309] = {
%      2,    3,    5,    7,    11,   13,   17,   19,   23,   29,   31,   37,
%      41,   43,   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,
%      97,   101,  103,  107,  109,  113,  127,  131,  137,  139,  149,  151,
%      157,  163,  167,  173,  179,  181,  191,  193,  197,  199,  211,  223,
%      227,  229,  233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
%      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,  353,  359,
%      367,  373,  379,  383,  389,  397,  401,  409,  419,  421,  431,  433,
%      439,  443,  449,  457,  461,  463,  467,  479,  487,  491,  499,  503,
%      509,  521,  523,  541,  547,  557,  563,  569,  571,  577,  587,  593,
%      599,  601,  607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
%      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,  739,  743,
%      751,  757,  761,  769,  773,  787,  797,  809,  811,  821,  823,  827,
%      829,  839,  853,  857,  859,  863,  877,  881,  883,  887,  907,  911,
%      919,  929,  937,  941,  947,  953,  967,  971,  977,  983,  991,  997,
%      1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
%      1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163,
%      1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249,
%      1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321,
%      1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439,
%      1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
%      1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601,
%      1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693,
%      1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783,
%      1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877,
%      1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
%      1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039};

%  // Row 52 in gold-127
%  static constexpr int gold_code[127] = {
%      1,  -1, 1,  -1, 1,  -1, 1,  -1, -1, -1, -1, -1, -1, -1, -1, -1,
%      1,  1,  -1, 1,  1,  1,  -1, 1,  -1, -1, -1, 1,  -1, -1, -1, -1,
%      -1, 1,  -1, 1,  -1, -1, 1,  -1, -1, 1,  1,  1,  -1, -1, 1,  1,
%      -1, 1,  1,  1,  1,  1,  -1, 1,  -1, 1,  1,  -1, -1, -1, 1,  1,
%      1,  1,  -1, 1,  -1, -1, 1,  1,  1,  -1, 1,  1,  -1, -1, -1, 1,
%      -1, -1, 1,  -1, -1, -1, -1, -1, -1, -1, 1,  -1, -1, 1,  1,  -1,
%      -1, -1, -1, 1,  -1, 1,  1,  1,  1,  -1, 1,  1,  -1, -1, 1,  1,
%      -1, 1,  -1, -1, 1,  -1, -1, 1,  -1, -1, 1,  1,  -1, -1, -1};
%};


% std::vector<size_t> CommsLib::getDataSc(size_t fftSize, size_t DataScNum,
%                                         size_t PilotScOffset) {
%   create data_sc (data subcarriers)    
%   std::vector<size_t> data_sc;
%   if (fftSize == Consts::kFftSize_80211) {
%
%     if size == 64 (this is the default and the case that we follow)
%     // We follow 802.11 PHY format here
%
%     assign index 0 - index 0 + kNumDataSubcarriers -> [0, 48] -> all provided subcarriers in lts_data_ind      
%     data_sc.assign(Consts::lts_data_ind,
%                    Consts::lts_data_ind + Consts::kNumDataSubcarriers_80211);
%   } else {  // Allocate the center subcarriers as data