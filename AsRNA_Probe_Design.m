function [results_array] = AsRNA_Probe_Design(rna_seq, prox_five_prime, ...
    prox_three_prime, anchor_five_prime, anchor_three_prime,...
    shortest_anchor_len, longest_anchor_len, shortest_prox_len,...
    longest_prox_len, max_CB_len)


%this code has been updated for the targeting of a tertiary contact site
%from a proximal "accessible anchor site" that will only be proximal if the
%tertiary contact is occuring 
%2 asRNAs will be designed based on user input: anchor (also referred to as
%TLR will be in position 1 (upstream of CB); proximal (also TL) will be in 
% position 2 (loop between CB and RBS)

%rna_seq is the sequence of the RNA of interest, inputted in between ' '
%prox_five_prime and prox_three_prime correspond to the nucleotide positions
%of the start and end (from 5' to 3') of the proposed long-range contact, respectively
%anchor_five_prime and anchor_three_prime correspond to the nucleotide positions
%of the start and end (from 5' to 3') of the proposed accessible anchor, respectively
%user can determine length of target regions (shortest_anchor_len_length, 
%longest_anchor_len_length, etc). Sequence to anchor target region will be added on 5'
%end. Sequence to proximal target region will be added on 3' end

% % %This is the P4-P6 region of the gI intron we are using:
% % %GAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTT
% % 
% % %The TetraLoop-TetraLoop Recepter matrix coordinates are, TL: 48th nucleotide to 51st, TLR: 122nd to 125th.
% % %TLR = [48,51]
% % %TL = [122,125]

%nucleotide positions for TL and TLR
TL_loc = [prox_five_prime:prox_three_prime];
TLR_loc = [anchor_five_prime:anchor_three_prime]; 

results = [];
results_array = []; %Creates an empty array.

%modifying length (and identity) of CB/RBS site
%all sequences are written 5' to 3'
CB_seq = 'TCACCTCTTGGAT';
%minimum nucleotide positions for the CB in order to preserve RBS
mandatory_CB = [1:8];
RBS_seq = 'ATTAAAGAGGAGA';
% note that mandatory_RBS = [6:13];

%if user inputs a max_CB_len that will not work with the design, print
%error
if max_CB_len > length(CB_seq) | max_CB_len < length(mandatory_CB);
    fprintf('max_CB_len input not possible \n');
end

counter=0;
%For each combination of anchor/proximal/CB length:
%Note that RBS length is always = CB length
   for anchor_length = shortest_anchor_len:longest_anchor_len;
       for prox_len = shortest_prox_len:longest_prox_len;
           for CB_len = length(mandatory_CB):max_CB_len;
                results = struct;
                counter = counter+1;
                [TLR_TL, probe, complete_probe, RBS_coords]=asRNA_design(rna_seq, TL_loc, TLR_loc, anchor_length, prox_len, CB_seq, RBS_seq, CB_len);
                %add columns corresponding to features of interest
                results.anchor_length = anchor_length; 
                results.prox_len = prox_len;
                results.CB_len = CB_len; 
                %'probe' results are the sequence from the asRNA 1 to end
                %of RBS
                results.probe=probe;
                results.target_1 = TLR_TL{1};
                results.target_2 = TLR_TL{2};                
                results.RBS_start = RBS_coords(:,1);
                results.RBS_end = RBS_coords(:,2);
                %complete_probe includes cutsite upstream as well as nucleotides
                %between end of RBS into GFP coding region (previously used
                %in biophysical model (vazquez-anderson et al)
                results.iRS3=complete_probe;
                
                %make files for nupack analysis
                iRS3_only = fopen(['nupack_commands/',num2str(counter),'.iRS3.in'],'wt');
                asRNA1_only = fopen(['nupack_commands/',num2str(counter),'.asRNA1.in'],'wt');
                asRNA2_only = fopen(['nupack_commands/',num2str(counter),'.asRNA2.in'],'wt');
                both_asRNA = fopen(['nupack_commands/',num2str(counter),'.both.in'],'wt');
                fprintf(iRS3_only, ['1 \n', probe, '\n1']); 
                fprintf(asRNA1_only, ['2 \n', probe, '\n',TLR_TL{1}, '\n1 2']); 
                fprintf(asRNA2_only, ['2 \n', probe, '\n', TLR_TL{2}, '\n1 2']); 
                fprintf(both_asRNA, ['3 \n', probe, '\n', TLR_TL{1}, '\n', TLR_TL{2}, '\n1 3 2']);
                fclose(iRS3_only);
                fclose(asRNA1_only);
                fclose(asRNA2_only);
                fclose(both_asRNA);
            results_array = [results_array, results];
            end
       end
   end

%    counter=0;
%    for 1:length(results_array);
%        counter = counter+1;
%         
%         commands = fopen(strcat('nupack_commands/commands.txt', 'w');
%         fprintf(commands, 'pairs -T 37 -multi -material rna target2
%     fprintf(fileID,'%6s %12s\n','x','exp(x)');
%     fprintf(fileID,'%6.2f %12.8f\n',A);
%     fclose(fileID);
%    end
end