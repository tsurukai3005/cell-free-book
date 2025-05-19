function [gainOverNoisedB_AP,gainOverNoisedB_BS,R_AP,R_BS,pilotIndex,BSassignment,D,D_small,APpositions,UEpositions,BSpositions,distances] = generateSetup_with_cellular(L,K,N,M,nbrOfSetups,seed,ASD_varphi,ASD_theta)
% この関数は、セクション5.3および5.4.3に記載されたシミュレーション設定の実現を生成します。
%
% 入力:
% L                  = セルフリーおよびスモールセルシステム用のAPの数
% K                  = ネットワーク内のUEの数
% N                  = セルフリーおよびスモールセルシステム用のAPごとのアンテナ数
% M                  = セルラーマッシブMIMO用のAPごとのアンテナ数
% nbrOfSetups        = ランダムなUEおよびAPの位置でのセットアップ数
% seed               = 乱数生成器のシード番号
% ASD_varphi         = 方位角における局所散乱モデルの角度標準偏差（ラジアン）
% ASD_theta          = 仰角における局所散乱モデルの角度標準偏差（ラジアン）
%
% 出力:
% gainOverNoisedB_AP = セルフリーおよびスモールセル設定で、AP l と UE k の間のセットアップ n のチャネルゲインの行列（ノイズ分散で正規化）
% gainOverNoisedB_BS = セルラーマッシブMIMO設定で、AP l と UE k の間のセットアップ n のチャネルゲインの行列（ノイズ分散で正規化）
% R_AP               = セルフリーおよびスモールセル設定で、AP l と UE k の間のセットアップ n の空間相関行列（ノイズ分散で正規化）
% R_BS               = セルラーマッシブMIMO設定で、AP l と UE k の間のセットアップ n の空間相関行列（ノイズ分散で正規化）
% pilotIndex         = UEに割り当てられたパイロットを含む、セットアップ数に対する行列
% BSassignment       = セルラーマッシブMIMOシステムでUEに割り当てられたAPのインデックスを含む、セットアップ数に対する行列
% D                  = DCC行列で、セルフリーセットアップにおいてAP l が UE k をサービスする場合は(l,k,n)が1、そうでなければ0
% D_small            = DCC行列で、スモールセルセットアップにおいてAP l が UE k をサービスする場合は(l,k,n)が1、そうでなければ0
% APpositions        = セルフリーおよびスモールセルシステム用のAPの位置を含む長さLのベクトル（実部は水平位置、虚部は垂直位置）
% UEpositions        = APpositionsと同じ方法で測定された長さKのUE位置を含むベクトル
% distances          = セルフリーおよびスモールセル設定におけるAPとUE間の距離を含む、gainOverNoisedB_APと同じ次元の行列
%
% このMatlab関数は、以下の研究成果をシミュレーションするために開発されました:
%
% Ozlem Tugfe Demir, Emil Bjornson, Luca Sanguinetti (2021),
% "Foundations of User-Centric Cell-Free Massive MIMO",
% Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
% pp 162-472. DOI: 10.1561/2000000109
%
% これはバージョン1.0です（最終編集日: 2021-01-31）
%
% ライセンス: このコードはGPLv2ライセンスの下でライセンスされています。このコードを使用して研究を行い、出版物につながる場合は、上記のモノグラフを引用してください。

%% シミュレーション設定を定義

% シード番号が0以外で指定されている場合は設定
if (nargin>5)&&(seed>0)
    rng(seed)
end

% カバレッジエリアのサイズ（ラップアラウンド付きの正方形）
squareLength = 1000; % メートル

% 通信帯域幅
B = 20e6;

% ノイズフィギュア（dB）
noiseFigure = 7;

% ノイズパワーを計算（dBm）
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

% (5.42)でのモデルのパスロスパラメータ
alpha = 36.7;
constantTerm = -30.5;

% (5.43)でのシャドウフェージングの標準偏差
sigma_sf = 4;

% シャドウフェージングの装飾距離
decorr = 9;

% APとUEの高さ差（メートル）
distanceVertical = 10;

% アンテナ間隔（波長の数）
antennaSpacing = 1/2; % 半波長距離

% セルラーAP（BS）の数
nbrBSs = 4;

% セルごとのUEの数と同じパイロット数
tau_p = K/nbrBSs;

% グリッド上のセルラーBSの数（次元ごと）
nbrBSsPerDim = sqrt(nbrBSs);

% BS間の距離（垂直/水平方向）
interBSDistance = squareLength/nbrBSsPerDim;

% グリッド上にBSを配置
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

% セルフリーおよびスモールセルセットアップ用のAPのグリッド上のAPの数
nbrAPsPerDim = sqrt(L);

% AP間の距離（垂直/水平方向）
interAPDistance = squareLength/nbrAPsPerDim;

% グリッド上にAPを配置
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

% ラップアラウンドを使用してBSとAPの代替位置を計算
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrBSs 1]);
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

% 結果を保存する準備
gainOverNoisedB_AP = zeros(L,K,nbrOfSetups);
gainOverNoisedB_BS = zeros(nbrBSs,K,nbrOfSetups);
distances = zeros(L,K,nbrOfSetups);
R_BS = zeros(M,M,nbrBSs,K,nbrOfSetups);
R_AP = zeros(N,N,L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);
BSassignment = zeros(K,nbrOfSetups);
masterAPs = zeros(K,1);
D = zeros(L,K,nbrOfSetups);
D_small = zeros(L,K,nbrOfSetups);

%% 全てのセットアップを通じて
for n = 1:nbrOfSetups

    % セルラーBSごとのUEの数を保存するために準備
    nbrOfUEsPerBS = zeros(nbrBSs,1);
    % ネットワーク内の落とされたUEの数を初期化
    nbrOfUEs = 0;
    % UEの位置を計算するために準備
    UEpositions = zeros(K,1);

    % シャドウフェージングの相関行列を保存するために準備
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowAPrealizations = zeros(K,L);
    shadowBSrealizations = zeros(K,nbrBSs);

    % UEを追加
    while nbrOfUEs<K

        % エリア内のランダムなUE位置を生成
        UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;

        % APがUEより10 m高いと仮定して距離を計算
        [distanceAPstoUE,whichposAP] = min(abs(APpositionsWrapped - repmat(UEposition,size(APpositionsWrapped))),[],2);
        distances(:,nbrOfUEs+1,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);

        % これが最初のUEでない場合
        if nbrOfUEs>0

            % 新しい予定のUEと他の全てのUEとの距離を計算
            shortestDistances = zeros(nbrOfUEs,1);

            for i = 1:nbrOfUEs
                shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
            end

            % 前のUEのシャドウフェージングの実現がすでに生成された場合に、新しいシャドウフェージングの実現を得るために必要な条件付き平均と標準偏差を計算
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:nbrOfUEs,1:nbrOfUEs);
            meanvaluesAP = term1*shadowAPrealizations(1:nbrOfUEs,:);
            meanvaluesBS = term1*shadowBSrealizations(1:nbrOfUEs,:);

            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);

        else % これが最初のUEの場合

            % UEを追加してシャドウフェージングの相関値を保存開始
            meanvaluesAP = 0;
            meanvaluesBS = 0;
            stdvalue = sigma_sf;
            newcolumn = [];

        end

        % セルフリーおよびスモールセル設定用のシャドウフェージングの実現を生成
        shadowing = meanvaluesAP + stdvalue*randn(1,L);

        % チャネルゲインをノイズパワーで割った値（dB）を計算
        gainsAP = constantTerm - alpha*log10(distances(:,nbrOfUEs+1,n)) + shadowing' - noiseVariancedBm;

        % 落とされたUEのマスターAPを最も良いチャネル条件のAPで決定
        [~,mastertemp] = max(gainsAP);

        % 最初のtau_p UEに直交パイロットを割り当て
        if nbrOfUEs+1 <= tau_p

            pilotIndextemp = nbrOfUEs+1;

        else % 残りのUEにパイロットを割り当て

            % マスターAPから各パイロットに対する受信パワーを計算
            pilotinterference = zeros(tau_p,1);

            for t = 1:tau_p

                pilotinterference(t) = sum(db2pow(gainOverNoisedB_AP(mastertemp,pilotIndex(1:nbrOfUEs,n)==t,n)));

            end

            % 最も受信パワーが少ないパイロットを見つける
            [~,bestpilot] = min(pilotinterference);
            pilotIndextemp = bestpilot;

        end

        % セルラーセットアップ用のシャドウフェージングの実現を生成

        shadowingBS = meanvaluesBS + stdvalue*randn(1,nbrBSs);


        % UEから各BSへの距離を計算
        [distanceBSstoUE,whichposBS] = min(abs(BSpositionsWrapped - repmat(UEposition,size(BSpositionsWrapped))),[],2);
        distancesBS = sqrt(distanceVertical^2+distanceBSstoUE.^2);

        % チャネルゲインをノイズパワーで割った値（dB）を計算
        gainsBS = constantTerm - alpha*log10(distancesBS) + shadowingBS' - noiseVariancedBm;

        % セルラーセットアップでUEが接続したいBSを見つける
        [~,bestBS] = max(gainsBS);

        % そのBSがまだtau_p UEを持っていない場合
        if nbrOfUEsPerBS(bestBS)<tau_p

            if sum(BSassignment(pilotIndex(:,n)==pilotIndextemp,n)==bestBS)<1


                % 好ましいセルにUEを追加
                nbrOfUEsPerBS(bestBS) = nbrOfUEsPerBS(bestBS) + 1;

                % これまでに落とされたUEの数を計算し、UEインデックスを保存
                nbrOfUEs = sum(nbrOfUEsPerBS);
                k = nbrOfUEs;

                % チャネルゲインをノイズパワーで割った値を保存
                gainOverNoisedB_AP(:,k,n) = gainsAP;
                gainOverNoisedB_BS(:,k,n) = gainsBS;


                % シャドウフェージングの相関行列を更新し、実現を保存
                shadowCorrMatrix(1:k-1,k) = newcolumn;
                shadowCorrMatrix(k,1:k-1) = newcolumn';
                shadowAPrealizations(k,:) = shadowing;
                shadowBSrealizations(k,:) = shadowingBS;

                % UEをセルラーBSに割り当て
                BSassignment(k,n) = bestBS;

                % UEの位置を保存
                UEpositions(k) = UEposition;


                % UE k のマスターAPを最も良いチャネル条件のAPで決定
                D(mastertemp,k,n) = 1;
                masterAPs(k) = mastertemp;

                % パイロット割り当てを保存
                pilotIndex(k,n) = pilotIndextemp;
                % 全てのAPを通じて
                for l = 1:L

                    % UE k と AP l の間の名目上の角度を計算
                    angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichposAP(l)));
                    angletoUE_theta = asin(distanceVertical/distances(l,k,n));
                    % 局所散乱モデルで空間相関行列を生成
                    if nargin>6
                        R_AP(:,:,l,k,n) = db2pow(gainOverNoisedB_AP(l,k,n))*functionRlocalscattering(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
                    else
                        R_AP(:,:,l,k,n) = db2pow(gainOverNoisedB_AP(l,k,n))*eye(N);  % 角度標準偏差が指定されていない場合、i.i.d. fading を設定
                    end
                end

                % 全てのBSを通じて
                for l = 1:nbrBSs

                    % 新しいUE k と BS l の間の名目上の角度を計算
                    angletoUE_varphi = angle(UEpositions(k)-BSpositionsWrapped(l,whichposBS(l)));
                    angletoUE_theta = asin(distanceVertical/distancesBS(l));
                    % 局所散乱モデルを使用して空間相関行列を生成
                    if nargin>6
                        R_BS(:,:,l,k,n) = db2pow(gainOverNoisedB_BS(l,k,n))*functionRlocalscattering(M,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
                    else
                        R_BS(:,:,l,k,n) = db2pow(gainOverNoisedB_BS(l,k,n))*eye(M);
                    end

                end



            end
        end

    end


    % セルフリーセットアップで、各パイロットに対して最も強いチャネル条件を持つUEにサービスを提供するAP
    for l = 1:L

        for t = 1:tau_p

            pilotUEs = find(t==pilotIndex(:,n));
            [~,UEindex] = max(gainOverNoisedB_AP(l,pilotUEs,n));
            D(l,pilotUEs(UEindex),n) = 1;



        end

    end

    % スモールセルセットアップにおいて、各UEにサービスを提供するAPを決定
    for k=1:K

        tempmat = -inf*ones(L,1);
        tempmat(D(:,k,n)==1,1) = gainOverNoisedB_AP(D(:,k,n)==1,k,n);
        [~,servingAP] = max(tempmat);
        D_small(servingAP,k,n) = 1;

    end


end
