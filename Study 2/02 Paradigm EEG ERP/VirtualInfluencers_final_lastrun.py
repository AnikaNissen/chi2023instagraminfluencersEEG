#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v3.2.4),
    on May 30, 2022, at 14:27
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '3.2.4'
expName = 'VirtualInfluencersEEG'  # from the Builder filename that created this script
expInfo = {'participant': '', 'session': '001'}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='C:\\Users\\eLearn\\Desktop\\Virtual Influencers\\VirtualInfluencers_final_lastrun.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(
    size=[1920, 1080], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=True, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='height')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "LSL"
LSLClock = core.Clock()
from psychopy import parallel
parallel.setPortAddress(61432) #61432
EndIntro = keyboard.Keyboard()

# Initialize components for Routine "StartExperiment"
StartExperimentClock = core.Clock()
StartNachricht1 = visual.TextStim(win=win, name='StartNachricht1',
    text='Press START to begin with the experiment.',
    font='Arial',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
Start = keyboard.Keyboard()

# Initialize components for Routine "conditions1"
conditions1Clock = core.Clock()
WebsiteImg = visual.ImageStim(
    win=win,
    name='WebsiteImg', 
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(0.55, 1),
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=512, interpolate=True, depth=-1.0)
Frage = visual.TextStim(win=win, name='Frage',
    text='default text',
    font='Arial',
    pos=(0, 0.1), height=0.06, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-2.0);
RatingBlau = visual.RatingScale(win=win, name='RatingBlau', low=1,
high=5,
scale='(Totally disagree)''(Totally agree)',
labels=['Totally disagree','','','','Totally agree'],
singleClick=True,
marker='triangle',
markerStart=3,
noMouse=True,
respKeys=['c','v','b','n','m'])

# Initialize components for Routine "Fixation"
FixationClock = core.Clock()
jitter = np.arange(2, 3, .25)
shuffle(jitter)
FixationCross = visual.ShapeStim(
    win=win, name='FixationCross', vertices='cross',
    size=(0.02, 0.02),
    ori=0, pos=(0, 0),
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-1.0, interpolate=True)

# Initialize components for Routine "Ende"
EndeClock = core.Clock()
text = visual.TextStim(win=win, name='text',
    text='\nDas Experiment ist beendet.\nSagen Sie der Versuchsleitung Bescheid.',
    font='Arial',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "LSL"-------
# update component parameters for each repeat
EndIntro.keys = []
EndIntro.rt = []
# keep track of which components have finished
LSLComponents = [EndIntro]
for thisComponent in LSLComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
LSLClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1
continueRoutine = True

# -------Run Routine "LSL"-------
while continueRoutine:
    # get current time
    t = LSLClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=LSLClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *EndIntro* updates
    waitOnFlip = False
    if EndIntro.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        EndIntro.frameNStart = frameN  # exact frame index
        EndIntro.tStart = t  # local t and not account for scr refresh
        EndIntro.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(EndIntro, 'tStartRefresh')  # time at next scr refresh
        EndIntro.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(EndIntro.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(EndIntro.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if EndIntro.status == STARTED and not waitOnFlip:
        theseKeys = EndIntro.getKeys(keyList=['z'], waitRelease=False)
        if len(theseKeys):
            theseKeys = theseKeys[0]  # at least one key was pressed
            
            # check for quit:
            if "escape" == theseKeys:
                endExpNow = True
            EndIntro.keys.append(theseKeys.name)  # storing all keys
            EndIntro.rt.append(theseKeys.rt)
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in LSLComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "LSL"-------
for thisComponent in LSLComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if EndIntro.keys in ['', [], None]:  # No response was made
    EndIntro.keys = None
thisExp.addData('EndIntro.keys',EndIntro.keys)
if EndIntro.keys != None:  # we had a response
    thisExp.addData('EndIntro.rt', EndIntro.rt)
thisExp.addData('EndIntro.started', EndIntro.tStartRefresh)
thisExp.addData('EndIntro.stopped', EndIntro.tStopRefresh)
thisExp.nextEntry()
# the Routine "LSL" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# ------Prepare to start Routine "StartExperiment"-------
# update component parameters for each repeat
Start.keys = []
Start.rt = []
# keep track of which components have finished
StartExperimentComponents = [StartNachricht1, Start]
for thisComponent in StartExperimentComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
StartExperimentClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1
continueRoutine = True

# -------Run Routine "StartExperiment"-------
while continueRoutine:
    # get current time
    t = StartExperimentClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=StartExperimentClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *StartNachricht1* updates
    if StartNachricht1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        StartNachricht1.frameNStart = frameN  # exact frame index
        StartNachricht1.tStart = t  # local t and not account for scr refresh
        StartNachricht1.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(StartNachricht1, 'tStartRefresh')  # time at next scr refresh
        StartNachricht1.setAutoDraw(True)
    
    # *Start* updates
    waitOnFlip = False
    if Start.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        Start.frameNStart = frameN  # exact frame index
        Start.tStart = t  # local t and not account for scr refresh
        Start.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(Start, 'tStartRefresh')  # time at next scr refresh
        Start.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(Start.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(Start.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if Start.status == STARTED and not waitOnFlip:
        theseKeys = Start.getKeys(keyList=['return'], waitRelease=False)
        if len(theseKeys):
            theseKeys = theseKeys[0]  # at least one key was pressed
            
            # check for quit:
            if "escape" == theseKeys:
                endExpNow = True
            Start.keys = theseKeys.name  # just the last key pressed
            Start.rt = theseKeys.rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in StartExperimentComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "StartExperiment"-------
for thisComponent in StartExperimentComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('StartNachricht1.started', StartNachricht1.tStartRefresh)
thisExp.addData('StartNachricht1.stopped', StartNachricht1.tStopRefresh)
# check responses
if Start.keys in ['', [], None]:  # No response was made
    Start.keys = None
thisExp.addData('Start.keys',Start.keys)
if Start.keys != None:  # we had a response
    thisExp.addData('Start.rt', Start.rt)
thisExp.addData('Start.started', Start.tStartRefresh)
thisExp.addData('Start.stopped', Start.tStopRefresh)
thisExp.nextEntry()
# the Routine "StartExperiment" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
trials = data.TrialHandler(nReps=1, method='fullRandom', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('conditions.xlsx'),
    seed=None, name='trials')
thisExp.addLoop(trials)  # add the loop to the experiment
thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
if thisTrial != None:
    for paramName in thisTrial:
        exec('{} = thisTrial[paramName]'.format(paramName))

for thisTrial in trials:
    currentLoop = trials
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            exec('{} = thisTrial[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "conditions1"-------
    # update component parameters for each repeat
    parallel.setData(int(marker))
    WebsiteImg.setImage(image)
    Frage.setText(question)
    RatingBlau.reset()
    # keep track of which components have finished
    conditions1Components = [WebsiteImg, Frage, RatingBlau]
    for thisComponent in conditions1Components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    conditions1Clock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    continueRoutine = True
    
    # -------Run Routine "conditions1"-------
    while continueRoutine:
        # get current time
        t = conditions1Clock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=conditions1Clock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *WebsiteImg* updates
        if WebsiteImg.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
            # keep track of start time/frame for later
            WebsiteImg.frameNStart = frameN  # exact frame index
            WebsiteImg.tStart = t  # local t and not account for scr refresh
            WebsiteImg.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(WebsiteImg, 'tStartRefresh')  # time at next scr refresh
            WebsiteImg.setAutoDraw(True)
        if WebsiteImg.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > WebsiteImg.tStartRefresh + 2.0-frameTolerance:
                # keep track of stop time/frame for later
                WebsiteImg.tStop = t  # not accounting for scr refresh
                WebsiteImg.frameNStop = frameN  # exact frame index
                win.timeOnFlip(WebsiteImg, 'tStopRefresh')  # time at next scr refresh
                WebsiteImg.setAutoDraw(False)
        
        # *Frage* updates
        if Frage.status == NOT_STARTED and tThisFlip >= 2.0-frameTolerance:
            # keep track of start time/frame for later
            Frage.frameNStart = frameN  # exact frame index
            Frage.tStart = t  # local t and not account for scr refresh
            Frage.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Frage, 'tStartRefresh')  # time at next scr refresh
            Frage.setAutoDraw(True)
        # *RatingBlau* updates
        if RatingBlau.status == NOT_STARTED and t >= 2.0-frameTolerance:
            # keep track of start time/frame for later
            RatingBlau.frameNStart = frameN  # exact frame index
            RatingBlau.tStart = t  # local t and not account for scr refresh
            RatingBlau.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(RatingBlau, 'tStartRefresh')  # time at next scr refresh
            RatingBlau.setAutoDraw(True)
        continueRoutine &= RatingBlau.noResponse  # a response ends the trial
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in conditions1Components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "conditions1"-------
    for thisComponent in conditions1Components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    trials.addData('WebsiteImg.started', WebsiteImg.tStartRefresh)
    trials.addData('WebsiteImg.stopped', WebsiteImg.tStopRefresh)
    trials.addData('Frage.started', Frage.tStartRefresh)
    trials.addData('Frage.stopped', Frage.tStopRefresh)
    # store data for trials (TrialHandler)
    trials.addData('RatingBlau.response', RatingBlau.getRating())
    trials.addData('RatingBlau.rt', RatingBlau.getRT())
    trials.addData('RatingBlau.started', RatingBlau.tStart)
    trials.addData('RatingBlau.stopped', RatingBlau.tStop)
    # the Routine "conditions1" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "Fixation"-------
    # update component parameters for each repeat
    jitter = np.arange(3.5, 4.5, .25)
    shuffle(jitter)
    # keep track of which components have finished
    FixationComponents = [FixationCross]
    for thisComponent in FixationComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    FixationClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    continueRoutine = True
    
    # -------Run Routine "Fixation"-------
    while continueRoutine:
        # get current time
        t = FixationClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=FixationClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *FixationCross* updates
        if FixationCross.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            FixationCross.frameNStart = frameN  # exact frame index
            FixationCross.tStart = t  # local t and not account for scr refresh
            FixationCross.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(FixationCross, 'tStartRefresh')  # time at next scr refresh
            FixationCross.setAutoDraw(True)
        if FixationCross.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > FixationCross.tStartRefresh + jitter[0]-frameTolerance:
                # keep track of stop time/frame for later
                FixationCross.tStop = t  # not accounting for scr refresh
                FixationCross.frameNStop = frameN  # exact frame index
                win.timeOnFlip(FixationCross, 'tStopRefresh')  # time at next scr refresh
                FixationCross.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in FixationComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "Fixation"-------
    for thisComponent in FixationComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    trials.addData('FixationCross.started', FixationCross.tStartRefresh)
    trials.addData('FixationCross.stopped', FixationCross.tStopRefresh)
    # the Routine "Fixation" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# completed 1 repeats of 'trials'


# ------Prepare to start Routine "Ende"-------
# update component parameters for each repeat
# keep track of which components have finished
EndeComponents = [text]
for thisComponent in EndeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
EndeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1
continueRoutine = True

# -------Run Routine "Ende"-------
while continueRoutine:
    # get current time
    t = EndeClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=EndeClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text* updates
    if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text.frameNStart = frameN  # exact frame index
        text.tStart = t  # local t and not account for scr refresh
        text.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
        text.setAutoDraw(True)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in EndeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "Ende"-------
for thisComponent in EndeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('text.started', text.tStartRefresh)
thisExp.addData('text.stopped', text.tStopRefresh)
# the Routine "Ende" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
