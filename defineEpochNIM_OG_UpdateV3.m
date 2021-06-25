function [eps] = defineEpochNIM_OG_UpdateV3(nantype)


names={'-MultiEnvSwitch','TMbase','Adapt (\DeltaEMG_{on(-)})','SplitPos','TR base','NIMBase',...
   'Task_{Switch}','Post1-Adapt_{SS}','Post1_{Early}','Post1_{Late}','Post2_{Early} - Post1_{Late}','Post2_{Late}','NIM Base','WithinEnvSwitc (-\DeltaEMG_{on(+)})'};

eps=defineEpochs(names,...
                {'OG base','TM tied 1','Neg Short','Pos Short','TR base','TR base'...
                'Adaptation','Post 1','Post 1','Post 1','Post 2','Post 2','TR base','TM tied 1',},...
                [-40 -40 20 20 -40 -40 -40 20 20 -40 20 -40 -40 -40],...
                [0,0,1,1,0,0,0,1,1,0,1,0, 0, 0],...
                [5,5,0,0,5,5,5,0,0,5,0,5, 5,5],...
                nantype);
   