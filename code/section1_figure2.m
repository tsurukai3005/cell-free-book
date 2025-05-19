%This Matlab script can be used to reproduce Figure 1.2 in the monograph:
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

% このMatlabスクリプトは、モノグラフの図1.2を再現するために使用できます：
% "ユーザー中心のセル・フリーマッシブMIMOの基礎"

%Empty workspace and close figures
% ワークスペースをクリアして図を閉じる
close all;
clear;

%Length of the coverage area in meter (as a square)
% カバレッジエリアの長さ（正方形、メートル単位）
squareLength = 1000;

%Total number of APs in the area
% エリア内のAP（アクセスポイント）の総数
nbrAPs = 4;

%Range of number of antennas and antenna gains
% アンテナ数とアンテナゲインの範囲
numberOfAntennas = [1 64];  % アンテナ数：1本または64本の二つのケース
antennaGains = [8 1];      % 対応するアンテナゲイン

%APs per dimension
% 次元ごとのAP数（グリッド状に配置するため）
nbrAPsPerDim = sqrt(nbrAPs);  % √4 = 2, 2×2のグリッド配置に

%Set the cell edge SNR to 0 dB) in case (a), which becomes -9 dB (1/8) when
%not accounting for the antenn gain
% セルエッジのSNRを0dBに設定（ケース(a)）、アンテナゲインを考慮しない場合は-9dB（1/8）
SNRedge = 1/8;

%Number of UE locations in horizontal/vertical dimension
% 水平/垂直方向のUE（ユーザー機器）位置の数
Kdims = 20;

%Number of UEs per cell
% セルごとのUE数
K = Kdims*Kdims;  % 20×20 = 400ユーザー/セル

%Pathloss exponent
% パスロス指数
alpha = 4;

%Maximum spectral efficiency
% 最大スペクトル効率
maxSE = 8;

%Bandwidth in MHz
% 帯域幅（MHz）
bandwidth = 10;


%Distance between APs in horizontal/vertical dimension
% 水平/垂直方向のAP間の距離
interSiteDistance = squareLength/nbrAPsPerDim;  % 1000/2 = 500メートル

%Put out the APs on a square grid
% APを正方形グリッド上に配置
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);  % 複素数で位置を表現（実部：x座標、虚部：y座標）


%Compute alternative AP locations by using wrap around
% ラップアラウンドを使用して代替AP位置を計算（境界効果を減らすため）
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrAPs 1]);
% ラップアラウンドにより、シミュレーションエリアを仮想的に拡張して境界のエッジ効果を軽減


%Prepare to put out the UEs in the cells
% セル内にUEを配置する準備
UEpositions = zeros(K,nbrAPs);  % 各UEの位置を格納する行列
perBS = zeros(nbrAPs,1);       % 各AP（基地局）のUE数をカウント

%Prepare to compute SEs
% スペクトル効率（SE）計算の準備
SEs = zeros(K,nbrAPs,length(numberOfAntennas));  % 各UE、各AP、各アンテナ設定ごとのSE値を格納


%Go through all the cells
% すべてのセルを処理
for l = 1:nbrAPs
    
    %Put out K users in each cell
    % 各セルにK人のユーザーを配置
    while min(perBS(l))<K
        
        posTmp = repmat(linspace(0,1,Kdims),[Kdims 1]);  % 0から1の範囲で等間隔に配置
        posTmp2 = posTmp';
        
        posX = interSiteDistance*posTmp(:) + real(BSpositions(l)) - interSiteDistance/2;  % x座標の計算
        posY = interSiteDistance*posTmp2(:) + imag(BSpositions(l)) - interSiteDistance/2; % y座標の計算
        
        UEpositions(:,l) = posX  + 1i*posY;  % 複素数形式で位置を格納
        perBS(l) = K;                       % このセルのUE数をKに設定
        
    end
    
    %Compute the distance from the UEs to AP l
    % UEからAPまでの距離を計算
    distancesSquaredBSl = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(l,:),[K 1]) ),[],2);
    % 複素数の絶対値（abs関数）を使用して距離を計算し、ラップアラウンドも考慮
    
    %Compute the signal power for each UE
    % 各UEの信号電力を計算
    signalPower = ((interSiteDistance/2)./distancesSquaredBSl).^(alpha);
    % パスロスモデル：(基準距離/実距離)^α
    
    %Compute the interference power for each UE
    % 各UEの干渉電力を計算
    interferencePower = zeros(size(distancesSquaredBSl));
    
    for j = 1:nbrAPs
        
        %Compute the distance from the users to UE j
        % ユーザーからAP jまでの距離を計算
        distancesSquaredBSj = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(j,:),[K 1]) ),[],2);
        
        %Add interference from non-serving APs
        % サービング中でないAPからの干渉を追加
        if j~=l
            
            interferencePower = interferencePower + ((interSiteDistance/2)./distancesSquaredBSj).^(alpha);
            % 他のAPからの干渉を合計
            
        end
        
    end
    
    %Compute SE for different number of antennas, assuming i.i.d. Rayleigh
    %fading and perfect CSI
    % 独立同一分布のレイリーフェージングと完全なCSI（チャネル状態情報）を仮定して
    % 異なるアンテナ数に対するSEを計算
    for m = 1:length(numberOfAntennas)
        
        SEs(:,l,m) = log2(1+numberOfAntennas(m)*antennaGains(m)*signalPower./(antennaGains(m)*interferencePower+1/SNRedge));
        % シャノン容量式を使用してSEを計算: log2(1 + SINR)
        % SINR = (アンテナ数×アンテナゲイン×信号電力)/(アンテナゲイン×干渉電力+ノイズ電力)
        
        SEs(SEs(:,l,m)>maxSE,l,m) = maxSE;
        % SEが最大値を超える場合は最大値に制限
        
    end
    
end

% 注：このスクリプトでは外部の関数を呼び出していません。すべての計算はスクリプト内で実行されています。
% 他のスクリプト（section1_figure10.mなど）では、computeSINRs_MMSE.mなどの関数を使用しています。

%Plot the simulation results
% シミュレーション結果をプロット
for m = 1:length(numberOfAntennas)
    
    figure;
    hold on; box on;
    
    for l = 1:nbrAPs
        
        surf(reshape(real(UEpositions(:,l)),[Kdims Kdims]),reshape(imag(UEpositions(:,l)),[Kdims Kdims]),bandwidth*reshape(SEs(:,l,m),[Kdims Kdims]));
        % 3Dサーフェスプロットを使用：x軸とy軸はUEの位置、z軸はデータレート（帯域幅×SE）
        
    end
    
    xlim([0 1000]);
    ylim([0 1000]);
    zlim([0 bandwidth*maxSE]);
    xlabel('Position [m]','Interpreter','latex');  % x軸ラベル：位置 [m]
    ylabel('Position [m]','Interpreter','latex');  % y軸ラベル：位置 [m]
    zlabel('Data rate [Mbit/s/user]','Interpreter','latex');  % z軸ラベル：データレート [Mbit/s/ユーザー]
    clim([0,bandwidth*maxSE]);
    view(-52,28);  % 3Dビューの角度を設定
    set(gca,'fontsize',16);  % フォントサイズを設定
    
end