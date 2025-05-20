%This Matlab script can be used to reproduce Figure 1.5 in the monograph:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

% このMatlabスクリプトは、モノグラフの図1.5を再現するために使用できます:
% 「ユーザー中心のセル・フリーマッシブMIMOの基礎」

%Empty workspace and close figures
% ワークスペースをクリアして図を閉じる
close all;
clear;

%Length of the coverage area in meter (as a square)
% カバレッジエリアの長さ（正方形、メートル単位）
squareLength = 1000;

%Total number of APs in the area
% エリア内のAP（アクセスポイント）の総数
nbrAPs = 64; % セル・フリー設定では64のAPを使用

%Range of number of antennas
% アンテナ数の範囲
numberOfAntennas = 1; % 各APに1本のアンテナを配置

%APs per dimension
% 次元ごとのAP数
nbrAPsPerDim = sqrt(nbrAPs); % 8×8のグリッド配置（√64 = 8）

%Set the cell edge SNR to -9 dB, which corresponds to 0 dB with an 9 dBi
%antenna as considered in Figure 1.2(a)
% セルエッジのSNRを-9dBに設定（図1.2(a)で考慮した9dBiアンテナでは0dBに相当）
SNRedge = 1/8; % -9dBは線形スケールで1/8

%Number of UE locations in horizontal/vertical dimension
% 水平/垂直方向のUE（ユーザー機器）位置の数
Kdims = 20;

%Number of UEs per cell
% セルごとのUE数
K = Kdims*Kdims; % 20×20=400ユーザー/セル

%Pathloss exponent
% パスロス指数
alpha = 4; % 都市部環境を想定した標準的な値

%Maximum spectral efficiency
% 最大スペクトル効率
maxSE = 8; % bit/s/Hz

%Bandwidth in MHz
% 帯域幅（MHz）
bandwidth = 10;


%Distance between APs in horizontal/vertical dimension
% 水平/垂直方向のAP間の距離（セル・フリー設定）
interSiteDistance = squareLength/nbrAPsPerDim; % 1000/8 = 125メートル

% セルラー設定でのAP間距離
interSiteDistanceCellular = squareLength/2; % 1000/2 = 500メートル


%Put out the APs on a square grid
% APを正方形グリッド上に配置
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
% 複素数を使用して2次元の位置を表現（実部がx座標、虚部がy座標）


%Compute alternative AP locations by using wrap around
% ラップアラウンドを使用して代替AP位置を計算（境界効果を軽減するため）
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrAPs 1]);
% 元の領域を周囲に8つコピーして拡張（3×3の領域構成）


%Prepare to put out the UEs in the cells
% セル内にUEを配置する準備
UEpositions = zeros(K,nbrAPs); % 各UEの位置を格納する行列
perBS = zeros(nbrAPs,1); % 各APのUE数をカウント

%Prepare to compute SEs
% スペクトル効率（SE）計算の準備
SEs_cellular = zeros(K,nbrAPs,length(numberOfAntennas)); % セルラー設定のSE
SEs_cellfree = zeros(K,nbrAPs,length(numberOfAntennas)); % セル・フリー設定のSE


%Go through all the cells
% すべてのセルを処理
for l = 1:nbrAPs
    
    %Put out K users in each cell
    % 各セルにK人のユーザーを配置
    while min(perBS(l))<K
        
        posTmp = repmat(linspace(0,1,Kdims),[Kdims 1]); % 正規化座標の格子点作成
        posTmp2 = posTmp';
        
        posX = interSiteDistance*posTmp(:) + real(BSpositions(l)) - interSiteDistance/2; % x座標計算
        posY = interSiteDistance*posTmp2(:) + imag(BSpositions(l)) - interSiteDistance/2; % y座標計算
        
        UEpositions(:,l) = posX  + 1i*posY; % 複素数形式で位置を格納
        perBS(l) = K; % このセルのUE数をKに設定
        
    end
    
    %Compute the distance from the UEs to AP l
    % UEからAP lまでの距離を計算
    distancesSquaredBSl = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(l,:),[K 1]) ),[],2);
    % ラップアラウンドを考慮して最小距離を選択
    
    %Compute the signal power for each UE
    % 各UEの信号電力を計算
    signalPower = ((interSiteDistanceCellular/2)./distancesSquaredBSl).^(alpha);
    % パスロスモデル: (基準距離/実際の距離)^α
    % 基準距離はセルラー設定の半径（interSiteDistanceCellular/2 = 250m）
    
    %Compute the interference power for each UE
    % 各UEの干渉電力を計算
    interferencePower = zeros(size(distancesSquaredBSl));
    
    for j = 1:nbrAPs
        
        %Compute the distance from the users to UE j
        % ユーザーからAP jまでの距離を計算
        distancesSquaredBSj = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(j,:),[K 1]) ),[],2);
        
        %Add interference from non-serving APs
        % サービス提供していないAPからの干渉を追加
        if j~=l
            
            interferencePower = interferencePower + ((interSiteDistanceCellular/2)./distancesSquaredBSj).^(alpha);
            % セルラー設定では干渉として扱う
            
        end
        
    end
    
    %Compute SE for different number of antennas, assuming i.i.d. Rayleigh
    %fading and perfect CSI
    % 独立同一分布のレイリーフェージングと完全なCSIを仮定して
    % 異なるアンテナ数に対するSEを計算
    for m = 1:length(numberOfAntennas)
        
        SEs_cellfree(:,l,m) = log2(1+SNRedge*numberOfAntennas(m)*(signalPower + interferencePower));
        % セル・フリー設定のSE計算: log₂(1 + SNR×アンテナ数×(信号電力+干渉電力))
        % セル・フリー設定では他セルの信号も協調的に利用するため、干渉が信号に変わる
        
        SEs_cellfree(SEs_cellfree(:,l,m)>maxSE,l,m) = maxSE;
        % 最大SEを制限
        
        SEs_cellular(:,l,m) = log2(1+numberOfAntennas(m)*signalPower./(interferencePower+1/SNRedge));
        % セルラー設定のSE計算: log₂(1 + アンテナ数×信号電力/(干渉電力+ノイズ))
        % 従来のセルラー方式では他セルの信号は干渉として扱われる
        
        SEs_cellular(SEs_cellular(:,l,m)>maxSE,l,m) = maxSE;
        % 最大SEを制限
        
    end
    
end


%Plot the simulation results
% シミュレーション結果をプロット
figure; % セル・フリー設定の結果を表示
hold on; box on;

for l = 1:nbrAPs
    
    surf(reshape(real(UEpositions(:,l)),[Kdims Kdims]),reshape(imag(UEpositions(:,l)),[Kdims Kdims]),bandwidth*reshape(SEs_cellfree(:,l,m),[Kdims Kdims]));
    % 3Dサーフェスプロット: x軸とy軸は位置、z軸はデータレート（帯域幅×SE）
    
end

xlim([0 1000]);
ylim([0 1000]);
zlim([0 bandwidth*maxSE]);
xlabel('Position [m]','Interpreter','latex'); % x軸ラベル
ylabel('Position [m]','Interpreter','latex'); % y軸ラベル
zlabel('Data rate [Mbit/s/user]','Interpreter','latex'); % z軸ラベル
clim([0,bandwidth*maxSE]); % カラースケール設定
view(-52,28); % 3Dビュー角度
set(gca,'fontsize',16); % フォントサイズ設定


figure; % セルラー設定の結果を表示
hold on; box on;

for l = 1:nbrAPs
    
    surf(reshape(real(UEpositions(:,l)),[Kdims Kdims]),reshape(imag(UEpositions(:,l)),[Kdims Kdims]),bandwidth*reshape(SEs_cellular(:,l,m),[Kdims Kdims]));
    % 3Dサーフェスプロット: x軸とy軸は位置、z軸はデータレート（帯域幅×SE）
    
end

xlim([0 1000]);
ylim([0 1000]);
zlim([0 bandwidth*maxSE]);
xlabel('Position [m]','Interpreter','latex'); % x軸ラベル
ylabel('Position [m]','Interpreter','latex'); % y軸ラベル
zlabel('Data rate [Mbit/s/user]','Interpreter','latex'); % z軸ラベル
clim([0,bandwidth*maxSE]); % カラースケール設定
view(-52,28); % 3Dビュー角度
set(gca,'fontsize',16); % フォントサイズ設定