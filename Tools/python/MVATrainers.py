import numpy
import pandas as pd
from DataGetter import DataGetter
from math import sqrt
from time import sleep
from glob import glob

def mainSKL(options):

  from sklearn.ensemble import RandomForestClassifier
  import xgboost as xgb
  import pickle

  print "PROCESSING TRAINING DATA"

  from taggerOptions import StandardVariables, getJetVarNames

  #get variables 
  globalVars, jetVars = StandardVariables(options.variables)
  allVars = globalVars + getJetVarNames(jetVars)

  # Import data
  dg = DataGetter(allVars)
  dataFiles = []
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_0_division_0_TTbarSingleLepT*_training_[01234].h5")
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_20_division_0_TTbarSingleLepT*_training_[01234].h5")
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_40_division_0_TTbarSingleLepT*_training_[01234].h5")
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_60_division_0_TTbarSingleLepT*_training_[01234].h5")
  print dataFiles
  trainData = dg.importData(samplesToRun = tuple(dataFiles), prescale=True, ptReweight=options.ptReweight)

  # Create random forest
  clf = xgb.XGBClassifier(base_score=0.5, colsample_bylevel=1, colsample_bytree=1,
                          gamma=0, learning_rate=0.001, max_delta_step=0, max_depth=6,
                          min_child_weight=0.1, missing=None, n_estimators=2000, nthread=28,
                          objective='binary:logistic', reg_alpha=0, reg_lambda=0.01,
                          scale_pos_weight=1, seed=0, silent=False, subsample=1 )
  #clf = RandomForestClassifier(n_estimators=500, max_depth=10, n_jobs = 28, verbose = True)

  print "TRAINING RF"
  
  # Train random forest 
  clf = clf.fit(trainData["data"], trainData["labels"][:,0], sample_weight=trainData["weights"][:,0])
  
  #Dump output from training
  fileObject = open(options.directory + "/" + "TrainingOutput.pkl",'wb')
  out = pickle.dump(clf, fileObject)
  fileObject.close()
      
  #output = clf.predict_proba(trainData["data"])[:,1]

def mainXGB(options):

  import xgboost as xgb
  from taggerOptions import StandardVariables, getJetVarNames

  print "PROCESSING TRAINING DATA"

  #get variables 
  globalVars, jetVars = StandardVariables(options.variables)
  allVars = globalVars + getJetVarNames(jetVars)

  # Import data
  dg = DataGetter(allVars)
  dataFiles = []
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_0_division_0_TTbarSingleLepT*_training_[01234].h5")
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_20_division_0_TTbarSingleLepT*_training_[01234].h5")
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_40_division_0_TTbarSingleLepT*_training_[01234].h5")
  dataFiles += glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_60_division_0_TTbarSingleLepT*_training_[01234].h5")
  trainData = dg.importData(samplesToRun = tuple(dataFiles), prescale=True, ptReweight=options.ptReweight)
  #= dg.importData(samplesToRun = tuple(glob(options.dataFilePath + "/trainingTuple_TTbarSingleLepT*_0_division_0_TTbarSingleLepT*_training_[0].h5")), prescale=True, ptReweight=options.ptReweight)

  print "TRAINING XGB"

  # Create xgboost classifier
  # Train random forest 
  xgData = xgb.DMatrix(trainData["data"], label=trainData["labels"][:,0])#, weight=trainData["weights"][:,0])
  param = {'max_depth':6, 'eta':0.03, 'objective':'binary:logistic', 'eval_metric':['error', 'auc', 'logloss'], 'nthread':28 }
  gbm = xgb.train(param, xgData, num_boost_round=2000)
  
  #Dump output from training
  gbm.save_model(options.directory + "/" + 'TrainingModel.xgb')

  #output = gbm.predict(xgData)


def mainTF(options):

  import tensorflow as tf
  from CreateModel import CreateModel
  from DataManager import DataManager
  from DataSet import DataSet

  print "PROCESSING TRAINING DATA"

  dgSig = DataGetter.DefinedVariables(options.netOp.vNames, signal = True)
  dgBg = DataGetter.DefinedVariables(options.netOp.vNames, background = True)

  print "Input Variables: ",len(dgSig.getList())

  # Import data
  print options.runOp.validationSamples
  validData = dgSig.importData(samplesToRun = tuple(options.runOp.validationSamples), ptReweight=options.runOp.ptReweight)
  validDataBG = dgBg.importData(samplesToRun = tuple(options.runOp.validationSamples), ptReweight=options.runOp.ptReweight)

  for key in validData:
    validData[key] = numpy.vstack([validData[key], validDataBG[key][:validData[key].shape[0]]])

  #get input/output sizes
  print validData["data"].shape
  nFeatures = validData["data"].shape[1]
  nLabels = validData["labels"].shape[1]
  nWeigts = validData["weights"].shape[1]

  #Training parameters
  l2Reg = options.runOp.l2Reg
  MiniBatchSize = options.runOp.minibatchSize
  nEpoch = options.runOp.nepoch
  ReportInterval = options.runOp.reportInterval
  validationCount = min(options.runOp.nValidationEvents, validData["data"].shape[0])

  #scale data inputs to mean 0, stddev 1
  categories = numpy.array(options.netOp.vCategories)
  mins = numpy.zeros(categories.shape, dtype=numpy.float32)
  ptps = numpy.zeros(categories.shape, dtype=numpy.float32)
  for i in xrange(categories.max()):
    selectedCategory = categories == i
    mins[selectedCategory] = validData["data"][:,selectedCategory].mean()
    ptps[selectedCategory] = validData["data"][:,selectedCategory].std()
  ptps[ptps < 1e-10] = 1.0

  ##Create data manager, this class controls how data is fed to the network for training
  #                 DataSet(fileGlob, xsec, Nevts, kFactor, sig, prescale, rescale)
  signalDataSets = [DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6p1/trainingTuple_*_division_0_TTbarSingleLepT_training_*.h5",      365.4,  61878989, 1.0, True,  1.0, 1.0, 5),
                    DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6p1/trainingTuple_*_division_0_TTbarSingleLepTbar_training_*.h5",   365.4,  61901450, 1.0, True,  1.0, 1.0, 5),]
  backgroundDataSets = [DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_TTbarSingleLepT_training_*.h5",    365.4,  61878989, 1.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_TTbarSingleLepTbar_training_*.h5", 365.4,  61901450, 1.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT100to200_training_*.h5",   27990000,  80684349, 0.0, False, 1.0, 1.0, 1), 
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT200to300_training_*.h5",   1712000 ,  57580393, 0.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT300to500_training_*.h5",   347700  ,  54537903, 0.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT500to700_training_*.h5",   32100   ,  62271343, 0.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT700to1000_training_*.h5",  6831    ,  45232316, 0.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT1000to1500_training_*.h5", 1207    ,  15127293, 0.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT1500to2000_training_*.h5", 119.9   ,  11826702, 0.0, False, 1.0, 1.0, 1),
                        DataSet("/cms/data/pastika/trainData_pt20_30_40_dRPi_tightMass_deepFlavor_v6/trainingTuple_*_division_0_QCD_HT2000toInf_training_*.h5",  25.24   ,   6039005, 0.0, False, 1.0, 1.0, 1),
                        ]

  dm = DataManager(options.netOp.vNames, nEpoch, nFeatures, nLabels, nWeigts, options.runOp.ptReweight, signalDataSets, backgroundDataSets)

  # Build the graph
  denseNetwork = [nFeatures]+options.netOp.denseLayers+[nLabels]
  convLayers = options.netOp.convLayers
  rnnNodes = options.netOp.rnnNodes
  rnnLayers = options.netOp.rnnLayers
  mlp = CreateModel(options, denseNetwork, convLayers, rnnNodes, rnnLayers, dm.inputDataQueue, MiniBatchSize, mins, 1.0/ptps)

  #summary writer
  summary_writer = tf.summary.FileWriter(options.runOp.directory + "log_graph", graph=tf.get_default_graph())

  print "TRAINING MLP"

  with tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=8) ) as sess:
    sess.run(tf.global_variables_initializer())

    #start queue runners
    dm.launchQueueThreads(sess)

    print "Reporting validation loss every %i batchces with %i events per batch for %i epochs"%(ReportInterval, MiniBatchSize, nEpoch)

    #preload the first data into staging area
    sess.run([mlp.stagingOp], feed_dict={mlp.reg: l2Reg, mlp.keep_prob:options.runOp.keepProb})

    i = 0
    try:
      while dm.continueTrainingLoop():
        #print sess.run(dm.inputDataQueue.dequeue_many(1000))[1]
        #run training operations 
        _, _, summary, _ = sess.run([mlp.stagingOp, mlp.train_step, mlp.merged_train_summary_op, mlp.batch_norm_ops], feed_dict={mlp.reg: l2Reg, mlp.keep_prob:options.runOp.keepProb, mlp.training: True})
        summary_writer.add_summary(summary, i)
        i += 1

        if i == 1 or not i % ReportInterval:
          #run validation operations 
          validation_loss, accuracy, summary_vl = sess.run([mlp.loss_ph, mlp.accuracy, mlp.merged_valid_summary_op], feed_dict={mlp.x_ph: validData["data"][:validationCount], mlp.y_ph_: validData["labels"][:validationCount], mlp.reg: l2Reg, mlp.wgt_ph: validData["weights"][:validationCount]})
          summary_writer.add_summary(summary_vl, i)
          print('Interval %d, validation accuracy %0.6f, validation loss %0.6f' % (i/ReportInterval, accuracy, validation_loss))

    except Exception, e:
      # Report exceptions to the coordinator.
      dm.requestStop(e)
    finally:
      # Terminate as usual. It is safe to call `coord.request_stop()` twice.
      dm.requestStop(e)
      dm.join()

    mlp.saveCheckpoint(sess, options.runOp.directory)
    mlp.saveModel(sess, options.runOp.directory)
