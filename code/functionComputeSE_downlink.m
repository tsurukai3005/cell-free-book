function [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
    SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
    Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR] ...
    = functionComputeSE_downlink(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot)
%Compute downlink SE for different transmit precoding schemes using the capacity
%bound in Theorem 6.1 for the centralized schemes and the capacity bound
%in Corollary 6.3 for the distributed schemes. Compute the genie-aided
%downlink SE from Corollary 6.6 for the centralized and the distributed operations.
%
% 異なる送信プリコーディングスキームを使用したダウンリンクSE（スペクトル効率）を計算します。
% 集中型スキームにはTheorem 6.1の容量境界、分散型スキームにはCorollary 6.3の
% 容量境界を使用します。また、集中型・分散型の両方について、Corollary 6.6に基づく
% genie-aided（完全CSI仮定）ダウンリンクSEも計算します。
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
% Hhat             = L*N × nbrOfRealizations × K の行列で、
%                    (:,n,k)はチャネル実現nにおけるUE kへの推定集合チャネル。
%
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
% H                = L*N × nbrOfRealizations × K の行列で、
%                    実際のチャネル実現を含みます。Hhatと同じ構造です。
%
%D                 = DCC matrix for cell-free setup with dimension L x K 
%                    where (l,k) is one if AP l serves UE k and zero otherwise
% D                = セル・フリー設定用のDCC（Dynamic Cooperation Clusters）行列、サイズはL×K。
%                   (l,k)は、AP lがUE kをサービスする場合は1、それ以外は0です。
%
%B                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel estimate 
%                    between AP l and UE k, normalized by noise variance
% B                = N×N×L×K の行列で、(:,:,l,k)はAP lとUE k間の
%                   チャネル推定の空間相関行列（ノイズ分散で正規化済み）。
%
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k,
%                    normalized by noise variance
% C                = N×N×L×K の行列で、(:,:,l,k)はAP lとUE k間の
%                   チャネル推定誤差の空間相関行列（ノイズ分散で正規化済み）。
%
%tau_c             = Length of coherence block
% tau_c            = コヒーレンスブロックの長さ
%
%tau_p             = Length of pilot sequences
% tau_p            = パイロットシーケンスの長さ
%
%nbrOfRealizations = Number of channel realizations
% nbrOfRealizations = チャネル実現の数（モンテカルロシミュレーションの反復回数）
%
%N                 = Number of antennas per AP
% N                = AP当たりのアンテナ数
%
%K                 = Number of UEs 
% K                = ユーザー機器（UE）の数
%
%L                 = Number of APs
% L                = アクセスポイント（AP）の数
%
%p                 = Uplink transmit power per UE (same for everyone)
% p                = UE当たりのアップリンク送信電力（全UEで同じ）
%
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k,
%                    normalized by noise
% R                = N×N×L×K の行列で、(:,:,l,k)はAP lとUE k間の
%                   空間相関行列（ノイズで正規化済み）。
%
%pilotIndex        = Vector containing the pilot assigned to each UE
% pilotIndex       = 各UEに割り当てられたパイロットを含むベクトル
%
%rho_dist          = Matrix with dimension L x K where (l,k) is the power
%                    allocated to UE k by AP l in the distributed downlink
%                    operation
% rho_dist         = L×K の行列で、(l,k)は分散ダウンリンク運用においてAP lがUE kに
%                   割り当てる電力を表します。
%
%gainOverNoisedB   = Matrix with dimension L x K where (l,k) is the channel
%                    gain (normalized by the noise variance) between AP l
%                    and UE k
% gainOverNoisedB  = L×K の行列で、(l,k)はAP lとUE k間のチャネルゲイン
%                   （ノイズ分散で正規化、dB単位）。
%
%rho_tot           = Maximum allowed transmit power for each AP 
% rho_tot          = 各APの最大許容送信電力
%
%OUTPUT:
%SE_MMSE           = SEs achieved with MMSE precoding in (6.16)
% SE_MMSE          = (6.16)のMMSEプリコーディングで達成されるSE値
%
%SE_P_MMSE         = SEs achieved with P-MMSE precoding in (6.17)
% SE_P_MMSE        = (6.17)のP-MMSEプリコーディングで達成されるSE値
%
%SE_P_RZF          = SEs achieved with P-RZF precoding in (6.18)
% SE_P_RZF         = (6.18)のP-RZFプリコーディングで達成されるSE値
%
%SE_L_MMSE         = SEs achieved with L-MMSE precoding in (6.25)
% SE_L_MMSE        = (6.25)のL-MMSEプリコーディングで達成されるSE値
%
%SE_LP_MMSE        = SEs achieved with LP-MMSE precoding in (6.33)
% SE_LP_MMSE       = (6.33)のLP-MMSEプリコーディングで達成されるSE値
%
%SE_MR             = SEs achieved with MR precoding in (6.26)
% SE_MR            = (6.26)のMR（Maximum Ratio）プリコーディングで達成されるSE値
%
%Gen_SE_P_MMSE     = Genie-aided SEs achieved with P-MMSE precoding in (6.17)
% Gen_SE_P_MMSE    = (6.17)のP-MMSEプリコーディングで達成されるgenie-aided SE値
%
%Gen_SE_P_RZF      = Genie-aided SEs achieved with P-RZF precoding in (6.18)
% Gen_SE_P_RZF     = (6.18)のP-RZFプリコーディングで達成されるgenie-aided SE値
%
%Gen_SE_LP_MMSE    = Genie-aided SEs achieved with LP-MMSE precoding in (6.33)
% Gen_SE_LP_MMSE   = (6.33)のLP-MMSEプリコーディングで達成されるgenie-aided SE値
%
%Gen_SE_MR         = Genie-aided SEs achieved with MR precoding in (6.26)
% Gen_SE_MR        = (6.26)のMRプリコーディングで達成されるgenie-aided SE値
%
%
%This Matlab function was developed to generate simulation results to:
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

% N×Nの単位行列を格納
eyeN = eye(N);

% ダウンリンクデータ送信のみを想定したプリログ因子を計算
prelogFactor = (1-tau_p/tau_c);

% シミュレーション結果の格納準備
Gen_SE_P_MMSE = zeros(K,1);
Gen_SE_P_RZF = zeros(K,1);
Gen_SE_LP_MMSE = zeros(K,1);
Gen_SE_MR = zeros(K,1);


% SE計算に使用される項の格納準備
signal_MR = zeros(K,1);
interf_MR = zeros(K,1);
cont_MR  = zeros(K,K);
scaling_MR = zeros(L,K);
interUserGains_MR = zeros(K,K,nbrOfRealizations);


signal_MMSE = zeros(K,1);
interf_MMSE = zeros(K,1);
scaling_MMSE = zeros(K,1);
portionScaling_MMSE = zeros(L,K);
interUserGains_MMSE = zeros(K,K,nbrOfRealizations);

signal_P_MMSE = zeros(K,1);
interf_P_MMSE = zeros(K,1);
scaling_P_MMSE = zeros(K,1);
portionScaling_PMMSE = zeros(L,K);
interUserGains_P_MMSE = zeros(K,K,nbrOfRealizations);

signal_P_RZF = zeros(K,1);
interf_P_RZF = zeros(K,1);
scaling_P_RZF = zeros(K,1);
portionScaling_PRZF = zeros(L,K);
interUserGains_P_RZF = zeros(K,K,nbrOfRealizations);

signal_L_MMSE = zeros(K,1);
interf_L_MMSE = zeros(K,1);
scaling_L_MMSE = zeros(L,K);
interUserGains_L_MMSE = zeros(K,K,nbrOfRealizations);

signal_LP_MMSE = zeros(K,1);
interf_LP_MMSE = zeros(K,1);
scaling_LP_MMSE = zeros(L,K);
interUserGains_LP_MMSE = zeros(K,K,nbrOfRealizations);


%% プリコーディングのスケーリング係数を計算
% MRプリコーディングの計算
for l = 1:L
    
    % AP lがサービスするUEを抽出
    servedUEs = find(D(l,:)==1);
    
    for ind = 1:length(servedUEs)
        
        % チャネル推定の空間相関行列を使用してスケーリング係数を計算
        scaling_MR(l,servedUEs(ind)) = trace(B(:,:,l,servedUEs(ind)));
        
    end
    
end

% すべてのチャネル実現をループ処理
for n = 1:nbrOfRealizations
    
    % L-MMSE、LP-MMSEプリコーディング
    for l = 1:L
        
        % APがサービスするUEを抽出
        servedUEs = find(D(l,:)==1);
        
        % AP lがサービスするUEの推定誤差共分散行列の合計を計算
        Cserved = sum(C(:,:,l,servedUEs),4);
        
        % L-MMSEとLP-MMSEプリコーディングを計算
        V_MR = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        V_L_MMSE = p*((p*(V_MR*V_MR')+p*sum(C(:,:,l,:),4)+eyeN)\V_MR(:,servedUEs));
        V_LP_MMSE = p*((p*(V_MR(:,servedUEs)*V_MR(:,servedUEs)')+p*Cserved+eyeN)\V_MR(:,servedUEs));
        
        % モンテカルロ法によりスケーリング係数を計算
        scaling_L_MMSE(l,servedUEs) = scaling_L_MMSE(l,servedUEs) + sum(abs(V_L_MMSE).^2,1)/nbrOfRealizations;
        scaling_LP_MMSE(l,servedUEs) = scaling_LP_MMSE(l,servedUEs) + sum(abs(V_LP_MMSE).^2,1)/nbrOfRealizations;
        
    end
    
    % MMSE、P-MMSE、P-RZFプリコーディング
    
    
    % すべてのUEをループ処理
    for k = 1:K
        
        % サービングAPの集合を決定
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);
        
        % UE kと部分的に同じAPセットによってサービスされるUEを判定
        % つまり、(5.15)の集合
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        
        % UE kのサービスに関わるAPのチャネル実現と推定誤差相関行列を抽出
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        
        % MMSE、P-MMSE、P-RZFプリコーディングを計算
        V_MMSE = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        V_P_MMSE = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        V_P_RZF = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));
        
        % モンテカルロ法によりスケーリング係数を計算
        scaling_MMSE(k) = scaling_MMSE(k) + sum(abs(V_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_MMSE(k) = scaling_P_MMSE(k) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_RZF(k) = scaling_P_RZF(k) + sum(abs(V_P_RZF).^2,1)/nbrOfRealizations;
        
        % すべてのサービングAPを処理
        for l=1:La
            
            % 集中型プリコーディングベクトルの部分を抽出
            V_MMSE2 = V_MMSE((l-1)*N+1:l*N,:);
            V_P_MMSE2 = V_P_MMSE((l-1)*N+1:l*N,:);
            V_P_RZF2 = V_P_RZF((l-1)*N+1:l*N,:);
            
            % Theorem 6.1の信号・干渉項の期待値内部の項の実現値を計算
            
            portionScaling_MMSE(servingAPs(l),k) = portionScaling_MMSE(servingAPs(l),k) ...
                + sum(abs(V_MMSE2).^2,1)/nbrOfRealizations;
            
            portionScaling_PMMSE(servingAPs(l),k) = portionScaling_PMMSE(servingAPs(l),k) ...
                + sum(abs(V_P_MMSE2).^2,1)/nbrOfRealizations;
            
            portionScaling_PRZF(servingAPs(l),k) = portionScaling_PRZF(servingAPs(l),k) ...
                + sum(abs(V_P_RZF2).^2,1)/nbrOfRealizations;
        end
    end
    
end

% 正規化集中型プリコーダの部分のノルム二乗を正規化
portionScaling_MMSE = portionScaling_MMSE./repmat(scaling_MMSE.',[L 1]);

portionScaling_PMMSE = portionScaling_PMMSE./repmat(scaling_P_MMSE.',[L 1]);

portionScaling_PRZF = portionScaling_PRZF./repmat(scaling_P_RZF.',[L 1]);


%% MR閉形式期待値を計算

% すべてのAPをループ処理
for l = 1:L
    
    % APがサービスするUEを抽出
    servedUEs = find(D(l,:)==1);
    
    % APがサービスするすべてのUEを処理
    for ind = 1:length(servedUEs)
        
        % UEインデックスを抽出
        k = servedUEs(ind);
        
        % (6.27)の所望信号項
        signal_MR(k) = signal_MR(k) + sqrt(rho_dist(l,k)*real(trace(B(:,:,l,k))));
        
        
        for i = 1:K
            
            
            % UE kからUE iへの非コヒーレント干渉（(6.28)の最初の項）
            interf_MR(i) = interf_MR(i) + rho_dist(l,k)*real(trace(B(:,:,l,k)*R(:,:,l,i)))/real(trace(B(:,:,l,k)));
            
            if pilotIndex(k) == pilotIndex(i)
                
                % UE kからUE iへのコヒーレント干渉（(6.28)の2番目の項）
                cont_MR(i,k) = cont_MR(i,k) + sqrt(rho_dist(l,k))*real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i)))/sqrt(real(trace(B(:,:,l,k))));
                
            end
            
        end
        
    end
    
end

% (7.43)のスケーラブルな集中ダウンリンク電力割り当てのパラメータ
upsilon = -0.5;
kappa = 0.5;

% (7.43)に基づく集中型プリコーディングの電力割り当て係数を計算
% 外部関数functionCentralizedPowerAllocationを呼び出し
rho_MMSE = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_MMSE,upsilon,kappa);

rho_PMMSE = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PMMSE,upsilon,kappa);

rho_PRZF = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PRZF,upsilon,kappa);

%% すべてのチャネル実現をループ処理
for n = 1:nbrOfRealizations
    
    
    
    % この実現のためのモンテカルロ結果を格納する行列
    interf_MMSE_n = zeros(K,K);
    interf_P_MMSE_n = zeros(K,K);
    interf_P_RZF_n = zeros(K,K);
    interf_L_MMSE_n = zeros(K,K);
    interf_LP_MMSE_n = zeros(K,K);
    
    % すべてのAPをループ処理
    for l = 1:L
        
        % すべてのUEからAP lへのチャネル実現を抽出
        Hallj = reshape(H((l-1)*N+1:l*N,n,:),[N K]);
        
        % すべてのUEからAP lへのチャネル推定を抽出
        Hhatallj = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        
        % AP lがサービスするUEを抽出
        servedUEs = find(D(l,:)==1);
        
        % AP lがサービスするUEの推定誤差共分散行列の合計を計算
        Cserved = sum(C(:,:,l,servedUEs),4);
        
        % MR結合（Maximum Ratio結合）を計算
        V_MR = Hhatallj(:,servedUEs);
        
        % L-MMSE結合を計算
        V_L_MMSE = p*((p*(Hhatallj*Hhatallj')+p*sum(C(:,:,l,:),4)+eyeN)\V_MR);
        
        % LP-MMSE結合を計算
        V_LP_MMSE = p*((p*(V_MR*V_MR')+p*Cserved+eyeN)\V_MR);
        
        
        % APがサービスするすべてのUEを処理
        for ind = 1:length(servedUEs)
            
            % UEインデックスを抽出
            k = servedUEs(ind);
            
            % MRプリコーディングを正規化
            w = V_MR(:,ind)*sqrt(rho_dist(l,k)/scaling_MR(l,k));
            
            % UEから他のUEに到達する信号のゲインを計算
            interUserGains_MR(:,k,n) = interUserGains_MR(:,k,n) + Hallj'*w;
            % L-MMSEプリコーディングを正規化
            w = V_L_MMSE(:,ind)*sqrt(rho_dist(l,k)/scaling_L_MMSE(l,k));
            
            % Corollary 6.3の信号・干渉項の期待値内部の項の実現値を計算
            signal_L_MMSE(k) = signal_L_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            interf_L_MMSE_n(:,k) = interf_L_MMSE_n(:,k) + Hallj'*w;
            
            % UEから他のUEに到達する信号のゲインを計算
            interUserGains_L_MMSE(:,k,n) = interUserGains_L_MMSE(:,k,n) + Hallj'*w;
            
            % LP-MMSEプリコーディングを正規化
            w = V_LP_MMSE(:,ind)*sqrt(rho_dist(l,k)/scaling_LP_MMSE(l,k));
            
            % Corollary 6.3の信号・干渉項の期待値内部の項の実現値を計算
            signal_LP_MMSE(k) = signal_LP_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            interf_LP_MMSE_n(:,k) = interf_LP_MMSE_n(:,k) + Hallj'*w;
            
            % UEから他のUEに到達する信号のゲインを計算
            interUserGains_LP_MMSE(:,k,n) = interUserGains_LP_MMSE(:,k,n) + Hallj'*w;
            
            
            
            
            
        end
        
    end
    
    
    
    % 集中型スキームの考慮
    
    
    % すべてのUEをループ処理
    for k = 1:K
        
        
        % サービングAPの集合を決定
        servingAPs = find(D(:,k)==1);
        
        La = length(servingAPs);
        
        % UE kと部分的に同じAPセットによってサービスされるUEを判定
        % つまり、(5.15)の集合
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        % UE kのサービスに関わるAPのチャネル実現と推定誤差相関行列を抽出
        Hallj_active = zeros(N*La,K);
        
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        % P-MMSEプリコーディングを計算
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));

        % 電力割り当てを適用
        w = w*sqrt(rho_PMMSE(k)/scaling_P_MMSE(k));
        
        
        % Theorem 6.1の信号・干渉項の期待値内部の項の実現値を計算
        signal_P_MMSE(k) = signal_P_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_MMSE_n(:,k) = interf_P_MMSE_n(:,k) + Hallj_active'*w;
        
        % UEから他のUEに到達する信号のゲインを計算
        interUserGains_P_MMSE(:,k,n) = interUserGains_P_MMSE(:,k,n) + Hallj_active'*w;
        
               
        % P-RZF結合を計算
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));
            
        % 電力割り当てを適用
        w = w*sqrt(rho_PRZF(k)/scaling_P_RZF(k));
        
        % Theorem 6.1の信号・干渉項の期待値内部の項の実現値を計算
        signal_P_RZF(k) = signal_P_RZF(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_RZF_n(:,k) = interf_P_RZF_n(:,k) + Hallj_active'*w;
        
        % UEから他のUEに到達する信号のゲインを計算
        interUserGains_P

        % UEから他のUEに到達する信号のゲインを計算
       interUserGains_P_RZF(:,k,n) = interUserGains_P_RZF(:,k,n) + Hallj_active'*w;
       
               
       % MMSE結合を計算
       w = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
           
       % 電力割り当てを適用
       w = w*sqrt(rho_MMSE(k)/scaling_MMSE(k));
       
       % Theorem 6.1の信号・干渉項の期待値内部の項の実現値を計算
       signal_MMSE(k) = signal_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
       interf_MMSE_n(:,k) = interf_MMSE_n(:,k) + Hallj_active'*w;
       
       % UEから他のUEに到達する信号のゲインを計算
       interUserGains_MMSE(:,k,n) = interUserGains_MMSE(:,k,n) + Hallj_active'*w;
       
       
       
   end
   
   % 1つの実現での干渉電力を計算
   interf_MMSE = interf_MMSE + sum(abs(interf_MMSE_n).^2,2)/nbrOfRealizations;
   interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;
   interf_P_RZF = interf_P_RZF + sum(abs(interf_P_RZF_n).^2,2)/nbrOfRealizations;
   
   interf_L_MMSE = interf_L_MMSE + sum(abs(interf_L_MMSE_n).^2,2)/nbrOfRealizations;
   interf_LP_MMSE = interf_LP_MMSE + sum(abs(interf_LP_MMSE_n).^2,2)/nbrOfRealizations;
   
   
   
end



%% SEを計算
% Corollary 6.4の閉形式表現を使ってCorollary 6.3のMRによるSEを計算
SE_MR = prelogFactor*real(log2(1+(abs(signal_MR).^2) ./ (interf_MR + sum(abs(cont_MR).^2,2) - abs(signal_MR).^2 + 1)));

% L-MMSEによるCorollary 6.3のSEを計算
SE_L_MMSE = prelogFactor*real(log2(1+(abs(signal_L_MMSE).^2) ./ (interf_L_MMSE - abs(signal_L_MMSE).^2 + 1)));

% LP-MMSEによるCorollary 6.3のSEを計算
SE_LP_MMSE = prelogFactor*real(log2(1+(abs(signal_LP_MMSE).^2) ./ (interf_LP_MMSE - abs(signal_LP_MMSE).^2 + 1)));

% MMSEによるTheorem 6.1のSEを計算
SE_MMSE = prelogFactor*real(log2(1+(abs(signal_MMSE).^2) ./ (interf_MMSE - abs(signal_MMSE).^2 + 1)));

% P-MMSEによるTheorem 6.1のSEを計算
SE_P_MMSE = prelogFactor*real(log2(1+(abs(signal_P_MMSE).^2) ./ (interf_P_MMSE - abs(signal_P_MMSE).^2 + 1)));

% P-RZFによるTheorem 6.1のSEを計算
SE_P_RZF = prelogFactor*real(log2(1+(abs(signal_P_RZF).^2) ./ (interf_P_RZF - abs(signal_P_RZF).^2 + 1)));




% Corollary 6.6のUEでの完全CSIを仮定したSE、すなわちgenie-aided SEを計算




% すべてのUEをループ処理
for k = 1:K
   
   % UEでの完全CSIを仮定したMRプリコーディングによるSEを計算
   
   Gen_SE_MR(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_MR(k,k,:)).^2 ./ ( sum(abs(interUserGains_MR(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
   
   
   % UEでの完全CSIを仮定したLP-MMSEプリコーディングによるSEを計算
   
   Gen_SE_LP_MMSE(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_LP_MMSE(k,k,:)).^2 ./ ( sum(abs(interUserGains_LP_MMSE(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
   
   
   % UEでの完全CSIを仮定したP-MMSEプリコーディングによるSEを計算
   
   Gen_SE_P_MMSE(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_P_MMSE(k,k,:)).^2 ./ ( sum(abs(interUserGains_P_MMSE(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
   
   % UEでの完全CSIを仮定したP-RZFプリコーディングによるSEを計算
   
   Gen_SE_P_RZF(k) = prelogFactor*mean(log2( 1 + abs(interUserGains_P_RZF(k,k,:)).^2 ./ ( sum(abs(interUserGains_P_RZF(k,[1:k-1 k+1:end],:)).^2,2) + 1) ),3);
   
end
% 未使用の大きな行列を削除
clear interUserGains_MR interUserGains_MMSE interUserGains_P_MMSE interUserGains_P_RZF interUserGains_L_MMSE interUserGains_LP_MMSE;