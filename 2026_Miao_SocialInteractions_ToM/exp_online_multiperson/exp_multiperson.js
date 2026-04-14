/***************************** 
 * Group3.4_Multiperson *
 *****************************/

import { core, data, sound, util, visual, hardware } from './lib/psychojs-2023.2.2.js';
const { PsychoJS } = core;
const { TrialHandler, MultiStairHandler } = data;
const { Scheduler } = util;
//some handy aliases as in the psychopy scripts;
const { abs, sin, cos, PI: pi, sqrt } = Math;
const { round } = util;


// store info about the experiment session:
let expName = 'group3.4_multiperson';  // from the Builder filename that created this script
let expInfo = {
    'participantID': '',
};

// Start code blocks for 'Before Experiment'
// init psychoJS:
const psychoJS = new PsychoJS({
  debug: true
});

// open window:
psychoJS.openWindow({
  fullscr: true,
  color: new util.Color([-1,-1,-1]),
  units: 'norm',
  waitBlanking: true,
  backgroundImage: '',
  backgroundFit: 'none',
});
// schedule the experiment:
psychoJS.schedule(psychoJS.gui.DlgFromDict({
  dictionary: expInfo,
  title: expName
}));

const flowScheduler = new Scheduler(psychoJS);
const dialogCancelScheduler = new Scheduler(psychoJS);
psychoJS.scheduleCondition(function() { return (psychoJS.gui.dialogComponent.button === 'OK'); }, flowScheduler, dialogCancelScheduler);

// flowScheduler gets run if the participants presses OK
flowScheduler.add(updateInfo); // add timeStamp
flowScheduler.add(experimentInit);
flowScheduler.add(consentRoutineBegin());
flowScheduler.add(consentRoutineEachFrame());
flowScheduler.add(consentRoutineEnd());
flowScheduler.add(opening1RoutineBegin());
flowScheduler.add(opening1RoutineEachFrame());
flowScheduler.add(opening1RoutineEnd());
flowScheduler.add(volumeAdjustRoutineBegin());
flowScheduler.add(volumeAdjustRoutineEachFrame());
flowScheduler.add(volumeAdjustRoutineEnd());
const instructionLoopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(instructionLoopLoopBegin(instructionLoopLoopScheduler));
flowScheduler.add(instructionLoopLoopScheduler);
flowScheduler.add(instructionLoopLoopEnd);


const audioLoopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(audioLoopLoopBegin(audioLoopLoopScheduler));
flowScheduler.add(audioLoopLoopScheduler);
flowScheduler.add(audioLoopLoopEnd);









flowScheduler.add(connectionRoutineBegin());
flowScheduler.add(connectionRoutineEachFrame());
flowScheduler.add(connectionRoutineEnd());
const textLoopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(textLoopLoopBegin(textLoopLoopScheduler));
flowScheduler.add(textLoopLoopScheduler);
flowScheduler.add(textLoopLoopEnd);









flowScheduler.add(endingRoutineBegin());
flowScheduler.add(endingRoutineEachFrame());
flowScheduler.add(endingRoutineEnd());
flowScheduler.add(quitPsychoJS, '', true);

// quit if user presses Cancel in dialog box:
dialogCancelScheduler.add(quitPsychoJS, '', false);

psychoJS.start({
  expName: expName,
  expInfo: expInfo,
  resources: [
    // resources:
    {'name': 'audioRunLoop_condition.xlsx', 'path': 'audioRunLoop_condition.xlsx'},
    {'name': 'audio/Narrative13Situation2.wav', 'path': 'audio/Narrative13Situation2.wav'},
    {'name': 'audio/Narrative7Situation1.wav', 'path': 'audio/Narrative7Situation1.wav'},
    {'name': 'audio/Narrative8Situation1.wav', 'path': 'audio/Narrative8Situation1.wav'},
    {'name': 'audio/Narrative5Situation1.wav', 'path': 'audio/Narrative5Situation1.wav'},
    {'name': 'audio/Narrative6Situation1.wav', 'path': 'audio/Narrative6Situation1.wav'},
    {'name': 'audio/Narrative7Situation2.wav', 'path': 'audio/Narrative7Situation2.wav'},
    {'name': 'audio/Narrative8Situation2.wav', 'path': 'audio/Narrative8Situation2.wav'},
    {'name': 'audio/Narrative5Situation2.wav', 'path': 'audio/Narrative5Situation2.wav'},
    {'name': 'audio/Narrative6Situation2.wav', 'path': 'audio/Narrative6Situation2.wav'},
    {'name': 'audio/Narrative7Situation3.wav', 'path': 'audio/Narrative7Situation3.wav'},
    {'name': 'audio/Narrative8Situation3.wav', 'path': 'audio/Narrative8Situation3.wav'},
    {'name': 'audio/Narrative5Situation3.wav', 'path': 'audio/Narrative5Situation3.wav'},
    {'name': 'audio/Narrative6Situation3.wav', 'path': 'audio/Narrative6Situation3.wav'},
    {'name': 'audio/Narrative7Situation4.wav', 'path': 'audio/Narrative7Situation4.wav'},
    {'name': 'audio/Narrative8Situation4.wav', 'path': 'audio/Narrative8Situation4.wav'},
    {'name': 'audio/Narrative5Situation4.wav', 'path': 'audio/Narrative5Situation4.wav'},
    {'name': 'audio/Narrative6Situation4.wav', 'path': 'audio/Narrative6Situation4.wav'},
    {'name': 'audio/Narrative7Situation5.wav', 'path': 'audio/Narrative7Situation5.wav'},
    {'name': 'audio/Narrative8Situation5.wav', 'path': 'audio/Narrative8Situation5.wav'},
    {'name': 'audio/Narrative5Situation5.wav', 'path': 'audio/Narrative5Situation5.wav'},
    {'name': 'audio/Narrative6Situation5.wav', 'path': 'audio/Narrative6Situation5.wav'},
    {'name': 'audio/Narrative7Situation6.wav', 'path': 'audio/Narrative7Situation6.wav'},
    {'name': 'audio/Narrative8Situation6.wav', 'path': 'audio/Narrative8Situation6.wav'},
    {'name': 'audio/Narrative5Situation6.wav', 'path': 'audio/Narrative5Situation6.wav'},
    {'name': 'audio/Narrative6Situation6.wav', 'path': 'audio/Narrative6Situation6.wav'},
    {'name': 'audio/Narrative7Situation7.wav', 'path': 'audio/Narrative7Situation7.wav'},
    {'name': 'audio/Narrative8Situation7.wav', 'path': 'audio/Narrative8Situation7.wav'},
    {'name': 'audio/Narrative5Situation7.wav', 'path': 'audio/Narrative5Situation7.wav'},
    {'name': 'audio/Narrative6Situation7.wav', 'path': 'audio/Narrative6Situation7.wav'},
    {'name': 'audio/Narrative7Situation8.wav', 'path': 'audio/Narrative7Situation8.wav'},
    {'name': 'audio/Narrative8Situation8.wav', 'path': 'audio/Narrative8Situation8.wav'},
    {'name': 'audio/Narrative5Situation8.wav', 'path': 'audio/Narrative5Situation8.wav'},
    {'name': 'audio/Narrative6Situation8.wav', 'path': 'audio/Narrative6Situation8.wav'},
    {'name': 'audio/Narrative7Situation9.wav', 'path': 'audio/Narrative7Situation9.wav'},
    {'name': 'audio/Narrative8Situation9.wav', 'path': 'audio/Narrative8Situation9.wav'},
    {'name': 'audio/Narrative5Situation9.wav', 'path': 'audio/Narrative5Situation9.wav'},
    {'name': 'audio/Narrative6Situation9.wav', 'path': 'audio/Narrative6Situation9.wav'},
    {'name': 'textRunLoop_condition.xlsx', 'path': 'textRunLoop_condition.xlsx'},
    {'name': 'image/consent.png', 'path': 'image/consent.png'},
    {'name': 'audio/Narrative13Situation1.wav', 'path': 'audio/Narrative13Situation1.wav'},
    {'name': 'default.png', 'path': 'https://pavlovia.org/assets/default/default.png'},
    {'name': 'image/procedure_continuous.png', 'path': 'image/procedure_continuous.png'}, 
    {'name': 'image/illustration_continuous.png', 'path': 'image/illustration_continuous.png'}, 
    {'name': 'image/illustration_continuous_circle.png', 'path': 'image/illustration_continuous_circle.png'}
  ]
});

psychoJS.experimentLogger.setLevel(core.Logger.ServerLevel.EXP);

// define a function that could generate random integer in [min, max)
function getRandomInt(min, max) {
  return Math.floor(Math.random() * (max - min)) + min;
}

var currentLoop;
var frameDur;
async function updateInfo() {
  currentLoop = psychoJS.experiment;  // right now there are no loops
  expInfo['date'] = util.MonotonicClock.getDateStr();  // add a simple timestamp
  expInfo['expName'] = expName;
  expInfo['psychopyVersion'] = '2023.2.2';
  expInfo['OS'] = window.navigator.platform;


  // store frame rate of monitor if we can measure it successfully
  expInfo['frameRate'] = psychoJS.window.getActualFrameRate();
  if (typeof expInfo['frameRate'] !== 'undefined')
    frameDur = 1.0 / Math.round(expInfo['frameRate']);
  else
    frameDur = 1.0 / 60.0; // couldn't get a reliable measure so guess

  // add info from the URL:
  util.addInfoFromUrl(expInfo);
  psychoJS.setRedirectUrls('https://app.prolific.com/submissions/complete?cc=C1JZXOLP', '');


  
  psychoJS.experiment.dataFileName = (("." + "/") + `data/${expInfo["participantID"]}_${expName}`);
  psychoJS.experiment.field_separator = '\t';


  return Scheduler.Event.NEXT;
}


var consentClock;
var consentImage;
var consentResponse;
var opening1Clock;
var parID;
var continuousQFull;
var nAudioRun;
var nTextRun;
var nInstruct;
var instructMes1;
var procedureImg1;
var scaleImg1;
var instructMes2;
var procedureImg2;
var scaleImg2;
var instructMes3;
var procedureImg3;
var scaleImg3;
var instructMes4;
var procedureImg4;
var scaleImg4;
var instructMes5;
var procedureImg5;
var scaleImg5;
var instructMes6;
var procedureImg6;
var scaleImg6;
var instructMes;
var procedureImg;
var scaleImg;
var instructText1;
var keyEndingOpening1;
var volumeAdjustClock;
var volumeText;
var endingAdjust;
var testingSound;
var opening2Clock;
var instructText;
var keyEndingOpening;
var mouse;
var illustrationProcedure;
var illustrationScale;
var readyForAudioClock;
var readyText_audio;
var fixationClock;
var fixISI;
var endPractice;
var audioPreClock;
var slider_decimals;
var slideSpeed;
var oldRating;
var slider_width;
var slider_height;
var slider_orientation;
var slider_granularity;
var tAudioStopped;
var audios;
var fixAudio;
var slider_shape;
var continuousSlider;
var questionText;
var questionExplanation;
var catchQuestionClock;
var catchQuestionText;
var catchQuestionSlider;
var reminderText_2;
var confirmResponse_2;
var intRatingClock;
var intQuestionText;
var intQuestionSlider;
var reminderText;
var confirmResponse;
var runIntervalBreakClock;
var breakText_audio;
var endingBreakKey;
var connectionClock;
var connectionInstruct;
var connectionEnding;
var readyForTextClock;
var readyText_text;
var textPreClock;
var texts;
var slider_shape_2;
var continuousSlider_2;
var questionText_2;
var questionExplanation_2;
var runBreakTextClock;
var breakText_text;
var endingBreakKey_2;
var endingClock;
var endingText;
var endingEnding;
var globalClock;
var routineTimer;
var frameTolerance;
frameTolerance = 0.001;

async function experimentInit() {
  // Initialize components for Routine "consent"
  consentClock = new util.Clock();
  consentImage = new visual.ImageStim({
    win : psychoJS.window,
    name : 'consentImage', units : undefined, 
    image : 'image/consent.png', mask : undefined,
    anchor : 'center',
    ori : 0.0, pos : [0, 0], size : [1.2, 1.5],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : 0.0 
  });
  consentResponse = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "opening1"
  opening1Clock = new util.Clock();
  // Run 'Begin Experiment' code from initializing1
  parID = Number.parseInt(expInfo["participantID"].split('').reverse().find(char => /\d/.test(char)), 10);
  continuousQFull = "(How much do you think there are multiple people in a scene at this moment?)";
  nAudioRun = 1;
  nTextRun = 1;
  nInstruct = 0;
  instructMes1 = "In this survey, you will listen to that voice reading out several narratives (short stories). After that, you will directly read additional narratives as their texts appear on the screen.\n\nYou will be invited to give your ratings about some contents in these stories. There are no right or wrong answers for these ratings.\n\nPress \"space\" to continue learning about the survey.";
  procedureImg1 = "image/procedure_continuous.png";
  scaleImg1 = "image/illustration_continuous.png";
  instructMes2 = "Every narrative is cut into **nine parts**. During each part, you will provide your ratings continuously at each moment, which are called **Continuous Ratings**. A slider will stay on the screen for the whole period of the narrative, as shown below (the left is when you listen to the narrative, and the right is when you read the narrative).\n\nPress \"space\" to continue, or press \"p\" to view the previous page.\n";
  procedureImg2 = "image/procedure_continuous.png";
  scaleImg2 = "image/illustration_continuous.png";
  instructMes3 = "To make continuous ratings, you would move your mouse on top of the slider, **without clicking**. A blue circle will indicate the position of your current rating, as shown below. Please change your ratings whenever you feel something different. The ratings will start as the narrative starts and stop a short time after the narrative ends.\n\nPress \"space\" to continue, or press \"p\" to view the previous page.";
  procedureImg3 = "image/procedure_continuous.png";
  scaleImg3 = "image/illustration_continuous_circle.png";
  instructMes4 = "The continuous ratings are about how much you think there are **multiple people** present in the current narrative scene at each moment. Specifically, you should give a high rating for \u201cmultiple people\u201d when there is more than one character and they are in the same physical location or scene. You should give a low rating for \u201cmultiple people\u201d when the current story is about only one character, or there is more than one character but they are not in the same physical location or scene.\n\nPress \"space\" to continue, or press \"p\" to view the previous page. \n";
  procedureImg4 = "image/procedure_continuous.png";
  scaleImg4 = "image/illustration_continuous.png";
  instructMes5 = "To give you some examples, the following sentences DO have multiple people:\n\u201cLinda told Amy about her situation.\u201d\n\u201cLucas saw Hank running towards him.\u201d\n\nThe following sentences DO NOT have multiple people:\n\u201cMark thought of his wife when he walked alone on the street.\u201d\n\u201cLucy heard a loud noise on her way home.\u201d\n\nPress \"space\" to continue, or press \"p\" to view the previous page.  \n";
  procedureImg5 = "image/procedure_continuous.png";
  scaleImg5 = "image/illustration_continuous.png";
  instructMes6 = "After each part of the story, you will answer another question about **on average**, how much you think there were multiple people in the last part of the story (as the orange bar above shows).\n\nTo help you get familiar with the task, at the beginning, one sample part of narratives will be used for practice. After that, all the following parts will be the formal survey.\n\nNow please press \u201cspace\u201d to begin practicing.";
  procedureImg6 = "image/procedure_continuous.png";
  scaleImg6 = "image/illustration_continuous.png";
  instructMes = [instructMes1, instructMes2, instructMes3, instructMes4, instructMes5, instructMes6];
  procedureImg = [procedureImg1, procedureImg2, procedureImg3, procedureImg4, procedureImg5, procedureImg6];
  scaleImg = [scaleImg1, scaleImg2, scaleImg3, scaleImg4, scaleImg5, scaleImg6];
  
  instructText1 = new visual.TextStim({
    win: psychoJS.window,
    name: 'instructText1',
    text: 'Thank you for participating in our survey!\n\nThis survey requires some audio to be played, so please ensure you have a working speaker or headphones connected to your computer.\n\nYou can adjust the volume to a comfortable level on the next page. \n\nPress “space” to continue.',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  keyEndingOpening1 = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "volumeAdjust"
  volumeAdjustClock = new util.Clock();
  volumeText = new visual.TextStim({
    win: psychoJS.window,
    name: 'volumeText',
    text: 'There should be an audio sample playing right now. \nPlease adjust the volume of your headphones or speaker to a comfortable value.\n\nIf your device is not playing any sounds, ensure it is fully connected to your computer. If it still does not work, you will not be able to complete the survey. \n\nPlease press "space" to continue, or press “Esc” if your audio is not working to quit the survey. ',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  endingAdjust = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  testingSound = new sound.Sound({
      win: psychoJS.window,
      value: 'audio/Narrative13Situation1.wav',
      secs: (- 1),
      });
  testingSound.setVolume(1.0);
  // Initialize components for Routine "opening2"
  opening2Clock = new util.Clock();
  instructText = new visual.TextStim({
    win: psychoJS.window,
    name: 'instructText',
    text: '',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  keyEndingOpening = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  mouse = new core.Mouse({
    win: psychoJS.window,
  });
  mouse.mouseClock = new util.Clock();
  illustrationProcedure = new visual.ImageStim({
    win : psychoJS.window,
    name : 'illustrationProcedure', units : undefined, 
    image : 'default.png', mask : undefined,
    anchor : 'center',
    ori : 0.0, pos : [0, 0.675], size : [1.6, 0.45],
    color : new util.Color([1,1,1]), opacity : 0.0,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -4.0 
  });
  illustrationScale = new visual.ImageStim({
    win : psychoJS.window,
    name : 'illustrationScale', units : undefined, 
    image : 'default.png', mask : undefined,
    anchor : 'center',
    ori : 0.0, pos : [0, (- 0.7)], size : [1.2, 0.5],
    color : new util.Color([1,1,1]), opacity : 0.0,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -5.0 
  });
  // Initialize components for Routine "readyForAudio"
  readyForAudioClock = new util.Clock();
  readyText_audio = new visual.TextStim({
    win: psychoJS.window,
    name: 'readyText_audio',
    text: 'Please get ready for listening to and rating the narratives.',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  // Initialize components for Routine "fixation"
  fixationClock = new util.Clock();
  // Run 'Begin Experiment' code from code_2
  /* Syntax Error: Fix Python code */
  fixISI = new visual.TextStim({
    win: psychoJS.window,
    name: 'fixISI',
    text: '',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.1,  wrapWidth: undefined, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  endPractice = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "audioPre"
  audioPreClock = new util.Clock();
  // Run 'Begin Experiment' code from code
  slider_decimals = 1;
  oldRating = 50;
  slider_width = 1;
  slider_height = 0.05;
  slider_orientation = 0;
  slider_granularity = 0.1;
  tAudioStopped = 0;

  if (frameDur <= 0.0833) {
    // to make the sampling rate relatively consistent across participants (83.3 ms per sample, for 60, 120, and 144 Hz frame rate)
    slideSpeed = Math.round(0.0833/frameDur);    
  } else {
    // if the frame rate is less than 12 Hz, then sample at every frame
    slideSpeed = 1
  }
  
  audios = new sound.Sound({
      win: psychoJS.window,
      value: 'A',
      secs: (- 1),
      });
  audios.setVolume(1.0);
  fixAudio = new visual.TextStim({
    win: psychoJS.window,
    name: 'fixAudio',
    text: '+',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.1,  wrapWidth: undefined, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -2.0 
  });
  
  slider_shape = new visual.Rect ({
    win: psychoJS.window, name: 'slider_shape', 
    width: [1.5, 0.15][0], height: [1.5, 0.15][1],
    ori: 0.0, pos: [0, (- 0.35)],
    anchor: 'center',
    lineWidth: 1.0, 
    colorSpace: 'rgb',
    lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: 0.0, depth: -3, interpolate: true,
  });
  
  continuousSlider = new visual.Slider({
    win: psychoJS.window, name: 'continuousSlider',
    startValue: undefined,
    size: [1.0, 0.1], pos: [0, (- 0.35)], ori: 0.0, units: psychoJS.window.units,
    labels: ["Not at all", "Partially", "Definitely yes"], fontSize: 0.07, ticks: [0, 50, 100],
    granularity: 0.0, style: ["RATING"],
    color: new util.Color([1.0, 1.0, 1.0]), markerColor: new util.Color([0.2157, 0.4118, 0.6706]), lineColor: new util.Color('White'), 
    opacity: undefined, fontFamily: 'Arial', bold: true, italic: false, depth: -4, 
    flip: false,
  });
  
  questionText = new visual.TextStim({
    win: psychoJS.window,
    name: 'questionText',
    text: 'Are there multiple people?',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0.5], height: 0.1,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -5.0 
  });

  
  questionExplanation = new visual.TextStim({
    win: psychoJS.window,
    name: 'questionExplanation',
    text: '',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0.35], height: 0.07,  wrapWidth: 1.8, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -7.0 
  });
  
  // Initialize components for Routine "catchQuestion"
  catchQuestionClock = new util.Clock();
  // Run 'Begin Experiment' code from end_this_question_2

  catchQuestionText = new visual.TextStim({
    win: psychoJS.window,
    name: 'catchQuestionText',
    text: '',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0.5], height: 0.1,  wrapWidth: 1.4, ori: 0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: 1,
    depth: -1.0 
  });
  
  catchQuestionSlider = new visual.Slider({
    win: psychoJS.window, name: 'catchQuestionSlider',
    startValue: undefined,
    size: [1.0, 0.1], pos: [0, (- 0.35)], ori: 0, units: psychoJS.window.units,
    labels: ["Not at all", "Partially", "Definitely yes"], fontSize: 0.07, ticks: [0, 50, 100],
    granularity: 0, style: ["RATING"],
    color: new util.Color('white'), markerColor: new util.Color('Red'), lineColor: new util.Color('White'), 
    opacity: 1, fontFamily: 'Arial', bold: true, italic: false, depth: -2, 
    flip: false,
  });
  
  reminderText_2 = new visual.TextStim({
    win: psychoJS.window,
    name: 'reminderText_2',
    text: '*Use your mouse to click the position where you think can best answer the question. Then press space to confirm and continue.',
    font: 'Calibri',
    units: undefined, 
    pos: [0, (- 0.8)], height: 0.07,  wrapWidth: 1.5, ori: 0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: 1,
    depth: -3.0 
  });
  
  confirmResponse_2 = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "intRating"
  intRatingClock = new util.Clock();
  intQuestionText = new visual.TextStim({
    win: psychoJS.window,
    name: 'intQuestionText',
    text: 'On average, how much do you think there were multiple people in a scene during the last story part?',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0.5], height: 0.1,  wrapWidth: 1.6, ori: 0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: 1,
    depth: -1.0 
  });
  
  intQuestionSlider = new visual.Slider({
    win: psychoJS.window, name: 'intQuestionSlider',
    startValue: undefined,
    size: [1.0, 0.1], pos: [0, (- 0.35)], ori: 0, units: psychoJS.window.units,
    labels: ["Not at all", "Partially", "Definitely yes"], fontSize: 0.07, ticks: [0, 50, 100],
    granularity: 0, style: ["RATING"],
    color: new util.Color('white'), markerColor: new util.Color('Red'), lineColor: new util.Color('White'), 
    opacity: 1, fontFamily: 'Arial', bold: true, italic: false, depth: -2, 
    flip: false,
  });
  
  reminderText = new visual.TextStim({
    win: psychoJS.window,
    name: 'reminderText',
    text: '*Use your mouse to click the position where you think can best answer the question. Then press space to confirm and continue.',
    font: 'Calibri',
    units: undefined, 
    pos: [0, (- 0.8)], height: 0.07,  wrapWidth: 1.5, ori: 0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: 1,
    depth: -3.0 
  });
  
  confirmResponse = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "runIntervalBreak"
  runIntervalBreakClock = new util.Clock();
  breakText_audio = new visual.TextStim({
    win: psychoJS.window,
    name: 'breakText_audio',
    text: 'Now you can break for a while (at most 1 minute).\n\nPlease press space to continue when you feel ready for the next part.',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  endingBreakKey = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "connection"
  connectionClock = new util.Clock();
  connectionInstruct = new visual.TextStim({
    win: psychoJS.window,
    name: 'connectionInstruct',
    text: 'For the next portion, you will read the narratives on the screen instead of hearing them. You will answer the same types of questions as you did for the previous narratives. Please read carefully, because the narratives will not stop playing and there is no opportunity to read the same part twice. As before, the first part of narrative is for you to practice.\n\nPress "space" to begin.',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  connectionEnding = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "readyForText"
  readyForTextClock = new util.Clock();
  readyText_text = new visual.TextStim({
    win: psychoJS.window,
    name: 'readyText_text',
    text: 'Please get ready for reading and rating the narratives.',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  // Initialize components for Routine "textPre"
  textPreClock = new util.Clock();
  texts = new visual.TextStim({
    win: psychoJS.window,
    name: 'texts',
    text: '',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0.1], height: 0.1,  wrapWidth: 1.9, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  slider_shape_2 = new visual.Rect ({
    win: psychoJS.window, name: 'slider_shape_2', 
    width: [1.5, 0.15][0], height: [1.5, 0.15][1],
    ori: 0.0, pos: [0, (- 0.35)],
    anchor: 'center',
    lineWidth: 1.0, 
    colorSpace: 'rgb',
    lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: 0.0, depth: -2, interpolate: true,
  });
  
  continuousSlider_2 = new visual.Slider({
    win: psychoJS.window, name: 'continuousSlider_2',
    startValue: undefined,
    size: [1.0, 0.1], pos: [0, (- 0.35)], ori: 0.0, units: psychoJS.window.units,
    labels: ["Not at all", "Partially", "Definitely yes"], fontSize: 0.07, ticks: [0, 50, 100],
    granularity: 0.0, style: ["RATING"],
    color: new util.Color([1.0, 1.0, 1.0]), markerColor: new util.Color([0.2157, 0.4118, 0.6706]), lineColor: new util.Color('White'), 
    opacity: undefined, fontFamily: 'Arial', bold: true, italic: false, depth: -3, 
    flip: false,
  });
  
  questionText_2 = new visual.TextStim({
    win: psychoJS.window,
    name: 'questionText_2',
    text: 'Are there multiple people?',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0.5], height: 0.1,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -4.0 
  });
  
  questionExplanation_2 = new visual.TextStim({
    win: psychoJS.window,
    name: 'questionExplanation_2',
    text: '',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0.35], height: 0.07,  wrapWidth: 1.8, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -6.0 
  });
  
  // Initialize components for Routine "runBreakText"
  runBreakTextClock = new util.Clock();
  breakText_text = new visual.TextStim({
    win: psychoJS.window,
    name: 'breakText_text',
    text: 'Now you can break for a while (at most 1 minute).\n\nPlease press space to continue when you feel ready for the next part.',
    font: 'Arial',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  endingBreakKey_2 = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Initialize components for Routine "ending"
  endingClock = new util.Clock();
  endingText = new visual.TextStim({
    win: psychoJS.window,
    name: 'endingText',
    text: 'Congratulations! You have completed the whole survey.\n\nThank you very much for your participation! You will receive your reimbursement for this within a week.\n\nPress "space" to end and quit.',
    font: 'Calibri',
    units: undefined, 
    pos: [0, 0], height: 0.08,  wrapWidth: 1.4, ori: 0.0,
    languageStyle: 'LTR',
    color: new util.Color('white'),  opacity: undefined,
    depth: -1.0 
  });
  
  endingEnding = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  // Create some handy timers
  globalClock = new util.Clock();  // to track the time since experiment started
  routineTimer = new util.CountdownTimer();  // to track time remaining of each (non-slip) routine
  
  return Scheduler.Event.NEXT;
}


var t;
var frameN;
var continueRoutine;
var _consentResponse_allKeys;
var consentComponents;
function consentRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'consent' ---
    t = 0;
    consentClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('consent.started', globalClock.getTime());
    consentResponse.keys = undefined;
    consentResponse.rt = undefined;
    _consentResponse_allKeys = [];
    // keep track of which components have finished
    consentComponents = [];
    consentComponents.push(consentImage);
    consentComponents.push(consentResponse);
    
    for (const thisComponent of consentComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function consentRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'consent' ---
    // get current time
    t = consentClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *consentImage* updates
    if (t >= 0.0 && consentImage.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      consentImage.tStart = t;  // (not accounting for frame time here)
      consentImage.frameNStart = frameN;  // exact frame index
      
      consentImage.setAutoDraw(true);
    }
    
    
    // *consentResponse* updates
    if (t >= 0.3 && consentResponse.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      consentResponse.tStart = t;  // (not accounting for frame time here)
      consentResponse.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { consentResponse.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { consentResponse.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { consentResponse.clearEvents(); });
    }
    
    if (consentResponse.status === PsychoJS.Status.STARTED) {
      let theseKeys = consentResponse.getKeys({keyList: ['y'], waitRelease: false});
      _consentResponse_allKeys = _consentResponse_allKeys.concat(theseKeys);
      if (_consentResponse_allKeys.length > 0) {
        consentResponse.keys = _consentResponse_allKeys[_consentResponse_allKeys.length - 1].name;  // just the last key pressed
        consentResponse.rt = _consentResponse_allKeys[_consentResponse_allKeys.length - 1].rt;
        consentResponse.duration = _consentResponse_allKeys[_consentResponse_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of consentComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function consentRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'consent' ---
    for (const thisComponent of consentComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('consent.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(consentResponse.corr, level);
    }
    psychoJS.experiment.addData('consentResponse.keys', consentResponse.keys);
    if (typeof consentResponse.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('consentResponse.rt', consentResponse.rt);
        psychoJS.experiment.addData('consentResponse.duration', consentResponse.duration);
        routineTimer.reset();
        }
    
    consentResponse.stop();
    // the Routine "consent" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var _keyEndingOpening1_allKeys;
var opening1Components;
function opening1RoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'opening1' ---
    t = 0;
    opening1Clock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('opening1.started', globalClock.getTime());
    // Run 'Begin Routine' code from initializing1
    instructText1.setAlignHoriz('left');
    
    keyEndingOpening1.keys = undefined;
    keyEndingOpening1.rt = undefined;
    _keyEndingOpening1_allKeys = [];
    // keep track of which components have finished
    opening1Components = [];
    opening1Components.push(instructText1);
    opening1Components.push(keyEndingOpening1);
    
    for (const thisComponent of opening1Components)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function opening1RoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'opening1' ---
    // get current time
    t = opening1Clock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *instructText1* updates
    if (t >= 0.0 && instructText1.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      instructText1.tStart = t;  // (not accounting for frame time here)
      instructText1.frameNStart = frameN;  // exact frame index
      
      instructText1.setAutoDraw(true);
    }
    
    
    // *keyEndingOpening1* updates
    if (t >= 0 && keyEndingOpening1.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      keyEndingOpening1.tStart = t;  // (not accounting for frame time here)
      keyEndingOpening1.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { keyEndingOpening1.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { keyEndingOpening1.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { keyEndingOpening1.clearEvents(); });
    }
    
    if (keyEndingOpening1.status === PsychoJS.Status.STARTED) {
      let theseKeys = keyEndingOpening1.getKeys({keyList: ['space'], waitRelease: false});
      _keyEndingOpening1_allKeys = _keyEndingOpening1_allKeys.concat(theseKeys);
      if (_keyEndingOpening1_allKeys.length > 0) {
        keyEndingOpening1.keys = _keyEndingOpening1_allKeys[_keyEndingOpening1_allKeys.length - 1].name;  // just the last key pressed
        keyEndingOpening1.rt = _keyEndingOpening1_allKeys[_keyEndingOpening1_allKeys.length - 1].rt;
        keyEndingOpening1.duration = _keyEndingOpening1_allKeys[_keyEndingOpening1_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of opening1Components)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function opening1RoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'opening1' ---
    for (const thisComponent of opening1Components) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('opening1.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(keyEndingOpening1.corr, level);
    }
    psychoJS.experiment.addData('keyEndingOpening1.keys', keyEndingOpening1.keys);
    if (typeof keyEndingOpening1.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('keyEndingOpening1.rt', keyEndingOpening1.rt);
        psychoJS.experiment.addData('keyEndingOpening1.duration', keyEndingOpening1.duration);
        routineTimer.reset();
        }
    
    keyEndingOpening1.stop();
    // the Routine "opening1" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var _endingAdjust_allKeys;
var volumeAdjustComponents;
function volumeAdjustRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'volumeAdjust' ---
    t = 0;
    volumeAdjustClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('volumeAdjust.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_12
    (endingAdjust.keys === "");
    volumeText.setAlignHoriz('left');
    
    endingAdjust.keys = undefined;
    endingAdjust.rt = undefined;
    _endingAdjust_allKeys = [];
    testingSound.setVolume(1.0);
    // keep track of which components have finished
    volumeAdjustComponents = [];
    volumeAdjustComponents.push(volumeText);
    volumeAdjustComponents.push(endingAdjust);
    volumeAdjustComponents.push(testingSound);
    
    for (const thisComponent of volumeAdjustComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function volumeAdjustRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'volumeAdjust' ---
    // get current time
    t = volumeAdjustClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // Run 'Each Frame' code from code_12
    if ((endingAdjust.keys === "escape")) {
        saveData({"thisExp": psychoJS.experiment});
        core.quit();
    }
    
    
    // *volumeText* updates
    if (t >= 0.0 && volumeText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      volumeText.tStart = t;  // (not accounting for frame time here)
      volumeText.frameNStart = frameN;  // exact frame index
      
      volumeText.setAutoDraw(true);
    }
    
    
    // *endingAdjust* updates
    if (t >= 0.3 && endingAdjust.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      endingAdjust.tStart = t;  // (not accounting for frame time here)
      endingAdjust.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { endingAdjust.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { endingAdjust.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { endingAdjust.clearEvents(); });
    }
    
    if (endingAdjust.status === PsychoJS.Status.STARTED) {
      let theseKeys = endingAdjust.getKeys({keyList: ['space', 'escape'], waitRelease: false});
      _endingAdjust_allKeys = _endingAdjust_allKeys.concat(theseKeys);
      if (_endingAdjust_allKeys.length > 0) {
        endingAdjust.keys = _endingAdjust_allKeys[_endingAdjust_allKeys.length - 1].name;  // just the last key pressed
        endingAdjust.rt = _endingAdjust_allKeys[_endingAdjust_allKeys.length - 1].rt;
        endingAdjust.duration = _endingAdjust_allKeys[_endingAdjust_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // start/stop testingSound
    if (t >= 0.5 && testingSound.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      testingSound.tStart = t;  // (not accounting for frame time here)
      testingSound.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ testingSound.play(); });  // screen flip
      testingSound.status = PsychoJS.Status.STARTED;
    }
    if (t >= (testingSound.getDuration() + testingSound.tStart)     && testingSound.status === PsychoJS.Status.STARTED) {
      testingSound.stop();  // stop the sound (if longer than duration)
      testingSound.status = PsychoJS.Status.FINISHED;
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of volumeAdjustComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function volumeAdjustRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'volumeAdjust' ---
    for (const thisComponent of volumeAdjustComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('volumeAdjust.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(endingAdjust.corr, level);
    }
    psychoJS.experiment.addData('endingAdjust.keys', endingAdjust.keys);
    if (typeof endingAdjust.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('endingAdjust.rt', endingAdjust.rt);
        psychoJS.experiment.addData('endingAdjust.duration', endingAdjust.duration);
        routineTimer.reset();
        }
    
    endingAdjust.stop();
    testingSound.stop();  // ensure sound has stopped at end of Routine
    // the Routine "volumeAdjust" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var instructionLoop;
function instructionLoopLoopBegin(instructionLoopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    instructionLoop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 999, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: undefined,
      seed: undefined, name: 'instructionLoop'
    });
    psychoJS.experiment.addLoop(instructionLoop); // add the loop to the experiment
    currentLoop = instructionLoop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    for (const thisInstructionLoop of instructionLoop) {
      snapshot = instructionLoop.getSnapshot();
      instructionLoopLoopScheduler.add(importConditions(snapshot));
      instructionLoopLoopScheduler.add(opening2RoutineBegin(snapshot));
      instructionLoopLoopScheduler.add(opening2RoutineEachFrame());
      instructionLoopLoopScheduler.add(opening2RoutineEnd(snapshot));
      instructionLoopLoopScheduler.add(instructionLoopLoopEndIteration(instructionLoopLoopScheduler, snapshot));
    }
    
    return Scheduler.Event.NEXT;
  }
}


async function instructionLoopLoopEnd() {
  // terminate loop
  psychoJS.experiment.removeLoop(instructionLoop);
  // update the current loop from the ExperimentHandler
  if (psychoJS.experiment._unfinishedLoops.length>0)
    currentLoop = psychoJS.experiment._unfinishedLoops.at(-1);
  else
    currentLoop = psychoJS.experiment;  // so we use addData from the experiment
  return Scheduler.Event.NEXT;
}


function instructionLoopLoopEndIteration(scheduler, snapshot) {
  // ------Prepare for next entry------
  return async function () {
    if (typeof snapshot !== 'undefined') {
      // ------Check if user ended loop early------
      if (snapshot.finished) {
        // Check for and save orphaned data
        if (psychoJS.experiment.isEntryEmpty()) {
          psychoJS.experiment.nextEntry(snapshot);
        }
        scheduler.stop();
      } else {
        psychoJS.experiment.nextEntry(snapshot);
      }
    return Scheduler.Event.NEXT;
    }
  };
}


var audioLoop;
function audioLoopLoopBegin(audioLoopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    audioLoop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 2, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: undefined,
      seed: undefined, name: 'audioLoop'
    });
    psychoJS.experiment.addLoop(audioLoop); // add the loop to the experiment
    currentLoop = audioLoop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    for (const thisAudioLoop of audioLoop) {
      snapshot = audioLoop.getSnapshot();
      audioLoopLoopScheduler.add(importConditions(snapshot));
      audioLoopLoopScheduler.add(readyForAudioRoutineBegin(snapshot));
      audioLoopLoopScheduler.add(readyForAudioRoutineEachFrame());
      audioLoopLoopScheduler.add(readyForAudioRoutineEnd(snapshot));
      const audioRunLoopLoopScheduler = new Scheduler(psychoJS);
      audioLoopLoopScheduler.add(audioRunLoopLoopBegin(audioRunLoopLoopScheduler, snapshot));
      audioLoopLoopScheduler.add(audioRunLoopLoopScheduler);
      audioLoopLoopScheduler.add(audioRunLoopLoopEnd);
      audioLoopLoopScheduler.add(runIntervalBreakRoutineBegin(snapshot));
      audioLoopLoopScheduler.add(runIntervalBreakRoutineEachFrame());
      audioLoopLoopScheduler.add(runIntervalBreakRoutineEnd(snapshot));
      audioLoopLoopScheduler.add(audioLoopLoopEndIteration(audioLoopLoopScheduler, snapshot));
    }
    
    return Scheduler.Event.NEXT;
  }
}


var audioRunLoop;
function audioRunLoopLoopBegin(audioRunLoopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    audioRunLoop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: 'audioRunLoop_condition.xlsx',
      seed: undefined, name: 'audioRunLoop'
    });
    psychoJS.experiment.addLoop(audioRunLoop); // add the loop to the experiment
    currentLoop = audioRunLoop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    for (const thisAudioRunLoop of audioRunLoop) {
      snapshot = audioRunLoop.getSnapshot();
      audioRunLoopLoopScheduler.add(importConditions(snapshot));
      audioRunLoopLoopScheduler.add(fixationRoutineBegin(snapshot));
      audioRunLoopLoopScheduler.add(fixationRoutineEachFrame());
      audioRunLoopLoopScheduler.add(fixationRoutineEnd(snapshot));
      audioRunLoopLoopScheduler.add(audioPreRoutineBegin(snapshot));
      audioRunLoopLoopScheduler.add(audioPreRoutineEachFrame());
      audioRunLoopLoopScheduler.add(audioPreRoutineEnd(snapshot));
      audioRunLoopLoopScheduler.add(intRatingRoutineBegin(snapshot));
      audioRunLoopLoopScheduler.add(intRatingRoutineEachFrame());
      audioRunLoopLoopScheduler.add(intRatingRoutineEnd(snapshot));
      audioRunLoopLoopScheduler.add(catchQuestionRoutineBegin(snapshot));
      audioRunLoopLoopScheduler.add(catchQuestionRoutineEachFrame());
      audioRunLoopLoopScheduler.add(catchQuestionRoutineEnd(snapshot));
      audioRunLoopLoopScheduler.add(audioRunLoopLoopEndIteration(audioRunLoopLoopScheduler, snapshot));
    }
    
    return Scheduler.Event.NEXT;
  }
}


async function audioRunLoopLoopEnd() {
  // terminate loop
  psychoJS.experiment.removeLoop(audioRunLoop);
  // update the current loop from the ExperimentHandler
  if (psychoJS.experiment._unfinishedLoops.length>0)
    currentLoop = psychoJS.experiment._unfinishedLoops.at(-1);
  else
    currentLoop = psychoJS.experiment;  // so we use addData from the experiment
  return Scheduler.Event.NEXT;
}


function audioRunLoopLoopEndIteration(scheduler, snapshot) {
  // ------Prepare for next entry------
  return async function () {
    if (typeof snapshot !== 'undefined') {
      // ------Check if user ended loop early------
      if (snapshot.finished) {
        // Check for and save orphaned data
        if (psychoJS.experiment.isEntryEmpty()) {
          psychoJS.experiment.nextEntry(snapshot);
        }
        scheduler.stop();
      } else {
        psychoJS.experiment.nextEntry(snapshot);
      }
    return Scheduler.Event.NEXT;
    }
  };
}


async function audioLoopLoopEnd() {
  // terminate loop
  psychoJS.experiment.removeLoop(audioLoop);
  // update the current loop from the ExperimentHandler
  if (psychoJS.experiment._unfinishedLoops.length>0)
    currentLoop = psychoJS.experiment._unfinishedLoops.at(-1);
  else
    currentLoop = psychoJS.experiment;  // so we use addData from the experiment
  return Scheduler.Event.NEXT;
}


function audioLoopLoopEndIteration(scheduler, snapshot) {
  // ------Prepare for next entry------
  return async function () {
    if (typeof snapshot !== 'undefined') {
      // ------Check if user ended loop early------
      if (snapshot.finished) {
        // Check for and save orphaned data
        if (psychoJS.experiment.isEntryEmpty()) {
          psychoJS.experiment.nextEntry(snapshot);
        }
        scheduler.stop();
      }
    return Scheduler.Event.NEXT;
    }
  };
}


var textLoop;
function textLoopLoopBegin(textLoopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    textLoop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 2, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: undefined,
      seed: undefined, name: 'textLoop'
    });
    psychoJS.experiment.addLoop(textLoop); // add the loop to the experiment
    currentLoop = textLoop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    for (const thisTextLoop of textLoop) {
      snapshot = textLoop.getSnapshot();
      textLoopLoopScheduler.add(importConditions(snapshot));
      textLoopLoopScheduler.add(readyForTextRoutineBegin(snapshot));
      textLoopLoopScheduler.add(readyForTextRoutineEachFrame());
      textLoopLoopScheduler.add(readyForTextRoutineEnd(snapshot));
      const textRunLoopLoopScheduler = new Scheduler(psychoJS);
      textLoopLoopScheduler.add(textRunLoopLoopBegin(textRunLoopLoopScheduler, snapshot));
      textLoopLoopScheduler.add(textRunLoopLoopScheduler);
      textLoopLoopScheduler.add(textRunLoopLoopEnd);
      textLoopLoopScheduler.add(runBreakTextRoutineBegin(snapshot));
      textLoopLoopScheduler.add(runBreakTextRoutineEachFrame());
      textLoopLoopScheduler.add(runBreakTextRoutineEnd(snapshot));
      textLoopLoopScheduler.add(textLoopLoopEndIteration(textLoopLoopScheduler, snapshot));
    }
    
    return Scheduler.Event.NEXT;
  }
}


var textRunLoop;
function textRunLoopLoopBegin(textRunLoopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    textRunLoop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: 'textRunLoop_condition.xlsx',
      seed: undefined, name: 'textRunLoop'
    });
    psychoJS.experiment.addLoop(textRunLoop); // add the loop to the experiment
    currentLoop = textRunLoop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    for (const thisTextRunLoop of textRunLoop) {
      snapshot = textRunLoop.getSnapshot();
      textRunLoopLoopScheduler.add(importConditions(snapshot));
      textRunLoopLoopScheduler.add(fixationRoutineBegin(snapshot));
      textRunLoopLoopScheduler.add(fixationRoutineEachFrame());
      textRunLoopLoopScheduler.add(fixationRoutineEnd(snapshot));
      textRunLoopLoopScheduler.add(textPreRoutineBegin(snapshot));
      textRunLoopLoopScheduler.add(textPreRoutineEachFrame());
      textRunLoopLoopScheduler.add(textPreRoutineEnd(snapshot));
      textRunLoopLoopScheduler.add(intRatingRoutineBegin(snapshot));
      textRunLoopLoopScheduler.add(intRatingRoutineEachFrame());
      textRunLoopLoopScheduler.add(intRatingRoutineEnd(snapshot));
      textRunLoopLoopScheduler.add(catchQuestionRoutineBegin(snapshot));
      textRunLoopLoopScheduler.add(catchQuestionRoutineEachFrame());
      textRunLoopLoopScheduler.add(catchQuestionRoutineEnd(snapshot));
      textRunLoopLoopScheduler.add(textRunLoopLoopEndIteration(textRunLoopLoopScheduler, snapshot));
    }
    
    return Scheduler.Event.NEXT;
  }
}


async function textRunLoopLoopEnd() {
  // terminate loop
  psychoJS.experiment.removeLoop(textRunLoop);
  // update the current loop from the ExperimentHandler
  if (psychoJS.experiment._unfinishedLoops.length>0)
    currentLoop = psychoJS.experiment._unfinishedLoops.at(-1);
  else
    currentLoop = psychoJS.experiment;  // so we use addData from the experiment
  return Scheduler.Event.NEXT;
}


function textRunLoopLoopEndIteration(scheduler, snapshot) {
  // ------Prepare for next entry------
  return async function () {
    if (typeof snapshot !== 'undefined') {
      // ------Check if user ended loop early------
      if (snapshot.finished) {
        // Check for and save orphaned data
        if (psychoJS.experiment.isEntryEmpty()) {
          psychoJS.experiment.nextEntry(snapshot);
        }
        scheduler.stop();
      } else {
        psychoJS.experiment.nextEntry(snapshot);
      }
    return Scheduler.Event.NEXT;
    }
  };
}


async function textLoopLoopEnd() {
  // terminate loop
  psychoJS.experiment.removeLoop(textLoop);
  // update the current loop from the ExperimentHandler
  if (psychoJS.experiment._unfinishedLoops.length>0)
    currentLoop = psychoJS.experiment._unfinishedLoops.at(-1);
  else
    currentLoop = psychoJS.experiment;  // so we use addData from the experiment
  return Scheduler.Event.NEXT;
}


function textLoopLoopEndIteration(scheduler, snapshot) {
  // ------Prepare for next entry------
  return async function () {
    if (typeof snapshot !== 'undefined') {
      // ------Check if user ended loop early------
      if (snapshot.finished) {
        // Check for and save orphaned data
        if (psychoJS.experiment.isEntryEmpty()) {
          psychoJS.experiment.nextEntry(snapshot);
        }
        scheduler.stop();
      }
    return Scheduler.Event.NEXT;
    }
  };
}


var instruction;
var scaleImage;
var procedureImage;
var _keyEndingOpening_allKeys;
var gotValidClick;
var opening2Components;
function opening2RoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'opening2' ---
    t = 0;
    opening2Clock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('opening2.started', globalClock.getTime());
    // Run 'Begin Routine' code from initializing2
    instructText.setAlignHoriz('left');
    instruction = instructMes[nInstruct];
    scaleImage = scaleImg[nInstruct];
    procedureImage = procedureImg[nInstruct];
    if ((nInstruct >= 1)) {
        illustrationProcedure.setOpacity(1);
    } else {
        illustrationProcedure.setOpacity(0);
    }
    if (((nInstruct >= 1) && (nInstruct <= 2))) {
        illustrationScale.setOpacity(1);
    } else {
        illustrationScale.setOpacity(0);
    }
    if (((nInstruct === 3) || (nInstruct === 4))) {
        instructText.setPos([0.0, (- 0.2)]);
    } else {
        instructText.setPos([0.0, (- 0.03)]);
    }
    
    instructText.setText(instruction);
    keyEndingOpening.keys = undefined;
    keyEndingOpening.rt = undefined;
    _keyEndingOpening_allKeys = [];
    // setup some python lists for storing info about the mouse
    // current position of the mouse:
    mouse.x = [];
    mouse.y = [];
    mouse.leftButton = [];
    mouse.midButton = [];
    mouse.rightButton = [];
    mouse.time = [];
    gotValidClick = false; // until a click is received
    illustrationProcedure.setImage(procedureImage);
    illustrationScale.setImage(scaleImage);
    // keep track of which components have finished
    opening2Components = [];
    opening2Components.push(instructText);
    opening2Components.push(keyEndingOpening);
    opening2Components.push(mouse);
    opening2Components.push(illustrationProcedure);
    opening2Components.push(illustrationScale);
    
    for (const thisComponent of opening2Components)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


var prevButtonState;
var _mouseButtons;
var _mouseXYs;
function opening2RoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'opening2' ---
    // get current time
    t = opening2Clock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *instructText* updates
    if (t >= 0.0 && instructText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      instructText.tStart = t;  // (not accounting for frame time here)
      instructText.frameNStart = frameN;  // exact frame index
      
      instructText.setAutoDraw(true);
    }
    
    
    // *keyEndingOpening* updates
    if (t >= 0 && keyEndingOpening.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      keyEndingOpening.tStart = t;  // (not accounting for frame time here)
      keyEndingOpening.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { keyEndingOpening.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { keyEndingOpening.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { keyEndingOpening.clearEvents(); });
    }
    
    if (keyEndingOpening.status === PsychoJS.Status.STARTED) {
      let theseKeys = keyEndingOpening.getKeys({keyList: ['space', 'p'], waitRelease: false});
      _keyEndingOpening_allKeys = _keyEndingOpening_allKeys.concat(theseKeys);
      if (_keyEndingOpening_allKeys.length > 0) {
        keyEndingOpening.keys = _keyEndingOpening_allKeys[_keyEndingOpening_allKeys.length - 1].name;  // just the last key pressed
        keyEndingOpening.rt = _keyEndingOpening_allKeys[_keyEndingOpening_allKeys.length - 1].rt;
        keyEndingOpening.duration = _keyEndingOpening_allKeys[_keyEndingOpening_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // *mouse* updates
    if (t >= 0.0 && mouse.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      mouse.tStart = t;  // (not accounting for frame time here)
      mouse.frameNStart = frameN;  // exact frame index
      
      mouse.status = PsychoJS.Status.STARTED;
      mouse.mouseClock.reset();
      prevButtonState = mouse.getPressed();  // if button is down already this ISN'T a new click
      }
    if (mouse.status === PsychoJS.Status.STARTED) {  // only update if started and not finished!
      _mouseButtons = mouse.getPressed();
      if (!_mouseButtons.every( (e,i,) => (e == prevButtonState[i]) )) { // button state changed?
        prevButtonState = _mouseButtons;
        if (_mouseButtons.reduce( (e, acc) => (e+acc) ) > 0) { // state changed to a new click
          _mouseXYs = mouse.getPos();
          mouse.x.push(_mouseXYs[0]);
          mouse.y.push(_mouseXYs[1]);
          mouse.leftButton.push(_mouseButtons[0]);
          mouse.midButton.push(_mouseButtons[1]);
          mouse.rightButton.push(_mouseButtons[2]);
          mouse.time.push(mouse.mouseClock.getTime());
        }
      }
    }
    
    // *illustrationProcedure* updates
    if (t >= 0.0 && illustrationProcedure.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      illustrationProcedure.tStart = t;  // (not accounting for frame time here)
      illustrationProcedure.frameNStart = frameN;  // exact frame index
      
      illustrationProcedure.setAutoDraw(true);
    }
    
    
    // *illustrationScale* updates
    if (t >= 0.0 && illustrationScale.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      illustrationScale.tStart = t;  // (not accounting for frame time here)
      illustrationScale.frameNStart = frameN;  // exact frame index
      
      illustrationScale.setAutoDraw(true);
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of opening2Components)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function opening2RoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'opening2' ---
    for (const thisComponent of opening2Components) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('opening2.stopped', globalClock.getTime());
    // Run 'End Routine' code from initializing2
    if ((keyEndingOpening.keys === "space")) {
        nInstruct = (nInstruct + 1);
    } else {
        if ((keyEndingOpening.keys === "p")) {
            nInstruct = (nInstruct - 1);
        }
    }
    if ((nInstruct < 0)) {
        nInstruct = 0;
    }
    if ((nInstruct === 6)) {
        instructionLoop.finished = true;
    }
    
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(keyEndingOpening.corr, level);
    }
    psychoJS.experiment.addData('keyEndingOpening.keys', keyEndingOpening.keys);
    if (typeof keyEndingOpening.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('keyEndingOpening.rt', keyEndingOpening.rt);
        psychoJS.experiment.addData('keyEndingOpening.duration', keyEndingOpening.duration);
        routineTimer.reset();
        }
    
    keyEndingOpening.stop();
    // store data for psychoJS.experiment (ExperimentHandler)
    psychoJS.experiment.addData('mouse.x', mouse.x);
    psychoJS.experiment.addData('mouse.y', mouse.y);
    psychoJS.experiment.addData('mouse.leftButton', mouse.leftButton);
    psychoJS.experiment.addData('mouse.midButton', mouse.midButton);
    psychoJS.experiment.addData('mouse.rightButton', mouse.rightButton);
    psychoJS.experiment.addData('mouse.time', mouse.time);
    
    // the Routine "opening2" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var nTrial;
var catchQ_presented;
var readyForAudioComponents;
function readyForAudioRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'readyForAudio' ---
    t = 0;
    readyForAudioClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    routineTimer.add(2.580000);
    // update component parameters for each repeat
    psychoJS.experiment.addData('readyForAudio.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_4
    nTrial = 0;
    catchQ_presented = false;
    
    // keep track of which components have finished
    readyForAudioComponents = [];
    readyForAudioComponents.push(readyText_audio);
    
    for (const thisComponent of readyForAudioComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


var frameRemains;
function readyForAudioRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'readyForAudio' ---
    // get current time
    t = readyForAudioClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *readyText_audio* updates
    if (t >= 0.0 && readyText_audio.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      readyText_audio.tStart = t;  // (not accounting for frame time here)
      readyText_audio.frameNStart = frameN;  // exact frame index
      
      readyText_audio.setAutoDraw(true);
    }
    
    frameRemains = 0.0 + 2.58 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (readyText_audio.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      readyText_audio.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of readyForAudioComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine && routineTimer.getTime() > 0) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function readyForAudioRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'readyForAudio' ---
    for (const thisComponent of readyForAudioComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('readyForAudio.stopped', globalClock.getTime());
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var isiText;
var continuousQExp;
var audioFile;
var fixDur;
var tFixationStart;
var _endPractice_allKeys;
var fixationComponents;
function fixationRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'fixation' ---
    t = 0;
    fixationClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('fixation.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_2
    if ((nAudioRun === 1)) {
        if ((nTrial === 0)) {
            isiText = "+";
            continuousQExp = continuousQFull;
        }
        if ((nTrial === 1)) {
            isiText = "You have completed the practice. Now please press space to go on to the formal survey.";
            continuousQExp = "";
        } else {
            if ((nTrial === 10)) {
                isiText = "New Story";
            } else {
                isiText = "+";
            }
        }
    } else {
        if ((nTrial === 1)) {
            isiText = "New Story";
        } else {
            if ((nTrial === 10)) {
                isiText = "New Story";
            } else {
                isiText = "+";
            }
        }
    }
    nTrial = (nTrial + 1);
    if ((((parID % 2) === 0) && (nAudioRun === 1))) {
        audioFile = conditions1_1;
        fixDur = dur1_1;
        psychoJS.experiment.addData('run', 1);
    } else {
        if ((((parID % 2) === 1) && (nAudioRun === 1))) {
            audioFile = conditions2_1;
            fixDur = dur2_1;
            psychoJS.experiment.addData('run', 1);
        } else {
            if ((((parID % 2) === 0) && (nAudioRun === 2))) {
                audioFile = conditions1_2;
                fixDur = dur1_2;
                psychoJS.experiment.addData('run', 2);
            } else {
                if ((((parID % 2) === 1) && (nAudioRun === 2))) {
                    audioFile = conditions2_2;
                    fixDur = dur2_2;
                    psychoJS.experiment.addData('run', 2);
                }
            }
        }
    }
    psychoJS.experiment.addData('audio_played', audioFile);
    tFixationStart = t;
    
    fixISI.setText(isiText);
    endPractice.keys = undefined;
    endPractice.rt = undefined;
    _endPractice_allKeys = [];
    // keep track of which components have finished
    fixationComponents = [];
    fixationComponents.push(fixISI);
    fixationComponents.push(endPractice);
    
    for (const thisComponent of fixationComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function fixationRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'fixation' ---
    // get current time
    t = fixationClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // Run 'Each Frame' code from code_2
    if (((nAudioRun === 2) && (nTrial === 1))) {
        continueRoutine = false;
    }
    if ((! ((nAudioRun === 1) && (nTrial === 2)))) {
        if ((t >= (tFixationStart + fixDur))) {
            continueRoutine = false;
        }
    }
    
    
    // *fixISI* updates
    if (t >= 0.0 && fixISI.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      fixISI.tStart = t;  // (not accounting for frame time here)
      fixISI.frameNStart = frameN;  // exact frame index
      
      fixISI.setAutoDraw(true);
    }
    
    
    // *endPractice* updates
    if (t >= 0.0 && endPractice.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      endPractice.tStart = t;  // (not accounting for frame time here)
      endPractice.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { endPractice.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { endPractice.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { endPractice.clearEvents(); });
    }
    
    if (endPractice.status === PsychoJS.Status.STARTED) {
      let theseKeys = endPractice.getKeys({keyList: ['space'], waitRelease: false});
      _endPractice_allKeys = _endPractice_allKeys.concat(theseKeys);
      if (_endPractice_allKeys.length > 0) {
        endPractice.keys = _endPractice_allKeys[_endPractice_allKeys.length - 1].name;  // just the last key pressed
        endPractice.rt = _endPractice_allKeys[_endPractice_allKeys.length - 1].rt;
        endPractice.duration = _endPractice_allKeys[_endPractice_allKeys.length - 1].duration;
        // a response ends the routine
        // only on practice trials
        if ((nAudioRun === 1) && (nTrial === 2)){
          continueRoutine = false;
        }
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of fixationComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function fixationRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'fixation' ---
    for (const thisComponent of fixationComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('fixation.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(endPractice.corr, level);
    }
    psychoJS.experiment.addData('endPractice.keys', endPractice.keys);
    if (typeof endPractice.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('endPractice.rt', endPractice.rt);
        psychoJS.experiment.addData('endPractice.duration', endPractice.duration);
        routineTimer.reset();
        }
    
    endPractice.stop();
    // the Routine "fixation" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var thisFrame;
var slider_data;
var mouseRec;
var audioStopped;
var audioPreComponents;
var mouseInSlider;

function audioPreRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'audioPre' ---
    t = 0;
    audioPreClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('audioPre.started', globalClock.getTime());
    // Run 'Begin Routine' code from code
    thisFrame = 0;
    slider_data = [];
    continuousSlider.markerPos = 50;
    mouseRec = mouse.getPos();
    audioStopped = false;
    mouseInSlider = false;
    
    audios.setValue(audioFile);
    audios.setVolume(1.0);
    continuousSlider.reset()
    questionExplanation.setText(continuousQExp);
    // keep track of which components have finished
    audioPreComponents = [];
    audioPreComponents.push(audios);
    audioPreComponents.push(fixAudio);
    audioPreComponents.push(slider_shape);
    audioPreComponents.push(continuousSlider);
    audioPreComponents.push(questionText);
    audioPreComponents.push(questionExplanation);
    
    for (const thisComponent of audioPreComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function audioPreRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'audioPre' ---
    // get current time
    t = audioPreClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // Run 'Each Frame' code from code
    
    
    if (((! slider_shape.contains(mouse)) && (! mouseInSlider))) {
      if (((thisFrame % slideSpeed) === 0)) {
          slider_data.push([999, Number.parseInt((t * 1000))]);
      }
    }
    if ((slider_shape.contains(mouse) && (! mouseInSlider))) {
        mouseInSlider = true;
    }
    if ((slider_shape.contains(mouse) && mouseInSlider)) {
        if ((mouse.getPos()[slider_orientation] !== mouseRec[slider_orientation])) {
            mouseRec = mouse.getPos();
            continuousSlider.markerPos = (((mouseRec[slider_orientation] / slider_width) * 100) + (100 / 2));
        }
        if (continuousSlider.markerPos) {
            if ((oldRating !== continuousSlider.markerPos)) {
                oldRating = continuousSlider.markerPos;
            }
        }
        if ((continuousSlider.markerPos && ((thisFrame % slideSpeed) === 0))) {
            slider_data.push([util.round(oldRating, slider_decimals), Number.parseInt((t * 1000))]);
        }
    }
    thisFrame = (thisFrame + 1);
    if (((audios.status === PsychoJS.Status.FINISHED) && (! audioStopped))) {
        tAudioStopped = t;
        audioStopped = true;
    }
    if (((audios.status === PsychoJS.Status.FINISHED) && (t >= ((0.1 + tAudioStopped) - frameTolerance)))) {
        continueRoutine = false;
    }
    
    // start/stop audios
    if (t >= 0.3 && audios.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      audios.tStart = t;  // (not accounting for frame time here)
      audios.frameNStart = frameN;  // exact frame index
      
      psychoJS.window.callOnFlip(function(){ audios.play(); });  // screen flip
      audios.status = PsychoJS.Status.STARTED;
    }
    if (t >= (audios.getDuration() + audios.tStart)     && audios.status === PsychoJS.Status.STARTED) {
      audios.stop();  // stop the sound (if longer than duration)
      audios.status = PsychoJS.Status.FINISHED;
    }
    
    // *fixAudio* updates
    if (t >= 0.0 && fixAudio.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      fixAudio.tStart = t;  // (not accounting for frame time here)
      fixAudio.frameNStart = frameN;  // exact frame index
      
      fixAudio.setAutoDraw(true);
    }
    
    
    // *slider_shape* updates
    if (t >= 0.0 && slider_shape.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      slider_shape.tStart = t;  // (not accounting for frame time here)
      slider_shape.frameNStart = frameN;  // exact frame index
      
      slider_shape.setAutoDraw(true);
    }
    
    
    // *continuousSlider* updates
    if (t >= 0.0 && continuousSlider.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      continuousSlider.tStart = t;  // (not accounting for frame time here)
      continuousSlider.frameNStart = frameN;  // exact frame index
      
      continuousSlider.setAutoDraw(true);
    }
    
    
    // *questionText* updates
    if (t >= 0.0 && questionText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      questionText.tStart = t;  // (not accounting for frame time here)
      questionText.frameNStart = frameN;  // exact frame index
      
      questionText.setAutoDraw(true);
    }
    
    
    // *questionExplanation* updates
    if (t >= 0.0 && questionExplanation.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      questionExplanation.tStart = t;  // (not accounting for frame time here)
      questionExplanation.frameNStart = frameN;  // exact frame index
      
      questionExplanation.setAutoDraw(true);
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of audioPreComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function audioPreRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'audioPre' ---
    for (const thisComponent of audioPreComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('audioPre.stopped', globalClock.getTime());
    // Run 'End Routine' code from code
    psychoJS.experiment.addData("continuousResponses_audio", slider_data);
    
    audios.stop();  // ensure sound has stopped at end of Routine
  
    // the Routine "audioPre" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var _confirmResponse_allKeys;
var intRatingComponents;
function intRatingRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'intRating' ---
    t = 0;
    intRatingClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('intRating.started', globalClock.getTime());
    // Run 'Begin Routine' code from end_this_question
    confirmResponse.keys = "";
    
    intQuestionSlider.reset()
    confirmResponse.keys = undefined;
    confirmResponse.rt = undefined;
    _confirmResponse_allKeys = [];
    // keep track of which components have finished
    intRatingComponents = [];
    intRatingComponents.push(intQuestionText);
    intRatingComponents.push(intQuestionSlider);
    intRatingComponents.push(reminderText);
    intRatingComponents.push(confirmResponse);
    
    for (const thisComponent of intRatingComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function intRatingRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'intRating' ---
    // get current time
    t = intRatingClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // Run 'Each Frame' code from end_this_question
    if (((nAudioRun === 2) && (nTrial === 1))) {
        continueRoutine = false;
    }
    if (((intQuestionSlider.history.length !== 0) && (confirmResponse.keys === "space"))) {
        continueRoutine = false;
    }
    
    
    // *intQuestionText* updates
    if (t >= 0.0 && intQuestionText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      intQuestionText.tStart = t;  // (not accounting for frame time here)
      intQuestionText.frameNStart = frameN;  // exact frame index
      
      intQuestionText.setAutoDraw(true);
    }
    
    
    // *intQuestionSlider* updates
    if (t >= 0.0 && intQuestionSlider.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      intQuestionSlider.tStart = t;  // (not accounting for frame time here)
      intQuestionSlider.frameNStart = frameN;  // exact frame index
      
      intQuestionSlider.setAutoDraw(true);
    }
    
    
    // *reminderText* updates
    if (t >= 0.0 && reminderText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      reminderText.tStart = t;  // (not accounting for frame time here)
      reminderText.frameNStart = frameN;  // exact frame index
      
      reminderText.setAutoDraw(true);
    }
    
    
    // *confirmResponse* updates
    if (t >= 0.5 && confirmResponse.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      confirmResponse.tStart = t;  // (not accounting for frame time here)
      confirmResponse.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { confirmResponse.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { confirmResponse.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { confirmResponse.clearEvents(); });
    }
    
    if (confirmResponse.status === PsychoJS.Status.STARTED) {
      let theseKeys = confirmResponse.getKeys({keyList: ['space'], waitRelease: false});
      _confirmResponse_allKeys = _confirmResponse_allKeys.concat(theseKeys);
      if (_confirmResponse_allKeys.length > 0) {
        confirmResponse.keys = _confirmResponse_allKeys[_confirmResponse_allKeys.length - 1].name;  // just the last key pressed
        confirmResponse.rt = _confirmResponse_allKeys[_confirmResponse_allKeys.length - 1].rt;
        confirmResponse.duration = _confirmResponse_allKeys[_confirmResponse_allKeys.length - 1].duration;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of intRatingComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function intRatingRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'intRating' ---
    for (const thisComponent of intRatingComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('intRating.stopped', globalClock.getTime());
    psychoJS.experiment.addData('intQuestionSlider.response', intQuestionSlider.getRating());
    psychoJS.experiment.addData('intQuestionSlider.rt', intQuestionSlider.getRT());
    psychoJS.experiment.addData('intQuestionSlider.history', intQuestionSlider.getHistory());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(confirmResponse.corr, level);
    }
    psychoJS.experiment.addData('confirmResponse.keys', confirmResponse.keys);
    if (typeof confirmResponse.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('confirmResponse.rt', confirmResponse.rt);
        psychoJS.experiment.addData('confirmResponse.duration', confirmResponse.duration);
        }
    
    confirmResponse.stop();
    // the Routine "intRating" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var random_number;
var random_number_2;
var catchQ;
var _confirmResponse_2_allKeys;
var catchQuestionComponents;
function catchQuestionRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'catchQuestion' ---
    t = 0;
    catchQuestionClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('catchQuestion.started', globalClock.getTime());
    // Run 'Begin Routine' code from end_this_question_2
    if (((nAudioRun === 2) && (nTrial === 1))) {
        continueRoutine = false;
    }
    confirmResponse_2.keys = "";
    random_number = getRandomInt(1, 16);
    if ((((random_number === 1) && (! catchQ_presented)) || ((nTrial === 19) && (! catchQ_presented)))) {
        catchQ_presented = true;
        if ((nTrial === 1)) {
            catchQ_presented = false;
            continueRoutine = false;
        }
    } else {
        continueRoutine = false;
    }
    random_number_2 = getRandomInt(1, 4);
    if ((random_number_2 === 1)) {
        catchQ = ("Please give a rating around 'Not at all' for this question.");
        psychoJS.experiment.addData("catch_label", "Not at all");
        psychoJS.experiment.addData("catch_value", 0);
    } else {
        if ((random_number_2 === 2)) {
            catchQ = ("Please give a rating around 'Partially' for this question.");
            psychoJS.experiment.addData("catch_label", 'Partially');
            psychoJS.experiment.addData("catch_value", 50);
        } else {
            if ((random_number_2 === 3)) {
                catchQ = ("Please give a rating around 'Definitely yes' for this question.");
                psychoJS.experiment.addData("catch_label", "Definitely yes");
                psychoJS.experiment.addData("catch_value", 100);
            }
        }
    }
    
    catchQuestionText.setText(catchQ);
    catchQuestionSlider.reset()
    confirmResponse_2.keys = undefined;
    confirmResponse_2.rt = undefined;
    _confirmResponse_2_allKeys = [];
    // keep track of which components have finished
    catchQuestionComponents = [];
    catchQuestionComponents.push(catchQuestionText);
    catchQuestionComponents.push(catchQuestionSlider);
    catchQuestionComponents.push(reminderText_2);
    catchQuestionComponents.push(confirmResponse_2);
    
    for (const thisComponent of catchQuestionComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function catchQuestionRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'catchQuestion' ---
    // get current time
    t = catchQuestionClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // Run 'Each Frame' code from end_this_question_2
    if (((catchQuestionSlider.history.length !== 0) && (confirmResponse_2.keys === "space"))) {
        continueRoutine = false;
    }
    
    
    // *catchQuestionText* updates
    if (t >= 0.0 && catchQuestionText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      catchQuestionText.tStart = t;  // (not accounting for frame time here)
      catchQuestionText.frameNStart = frameN;  // exact frame index
      
      catchQuestionText.setAutoDraw(true);
    }
    
    
    // *catchQuestionSlider* updates
    if (t >= 0.0 && catchQuestionSlider.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      catchQuestionSlider.tStart = t;  // (not accounting for frame time here)
      catchQuestionSlider.frameNStart = frameN;  // exact frame index
      
      catchQuestionSlider.setAutoDraw(true);
    }
    
    
    // *reminderText_2* updates
    if (t >= 0.0 && reminderText_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      reminderText_2.tStart = t;  // (not accounting for frame time here)
      reminderText_2.frameNStart = frameN;  // exact frame index
      
      reminderText_2.setAutoDraw(true);
    }
    
    
    // *confirmResponse_2* updates
    if (t >= 0.5 && confirmResponse_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      confirmResponse_2.tStart = t;  // (not accounting for frame time here)
      confirmResponse_2.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { confirmResponse_2.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { confirmResponse_2.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { confirmResponse_2.clearEvents(); });
    }
    
    if (confirmResponse_2.status === PsychoJS.Status.STARTED) {
      let theseKeys = confirmResponse_2.getKeys({keyList: ['space'], waitRelease: false});
      _confirmResponse_2_allKeys = _confirmResponse_2_allKeys.concat(theseKeys);
      if (_confirmResponse_2_allKeys.length > 0) {
        confirmResponse_2.keys = _confirmResponse_2_allKeys[_confirmResponse_2_allKeys.length - 1].name;  // just the last key pressed
        confirmResponse_2.rt = _confirmResponse_2_allKeys[_confirmResponse_2_allKeys.length - 1].rt;
        confirmResponse_2.duration = _confirmResponse_2_allKeys[_confirmResponse_2_allKeys.length - 1].duration;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of catchQuestionComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function catchQuestionRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'catchQuestion' ---
    for (const thisComponent of catchQuestionComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('catchQuestion.stopped', globalClock.getTime());
    psychoJS.experiment.addData('catchQuestionSlider.response', catchQuestionSlider.getRating());
    psychoJS.experiment.addData('catchQuestionSlider.rt', catchQuestionSlider.getRT());
    psychoJS.experiment.addData('catchQuestionSlider.history', catchQuestionSlider.getHistory());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(confirmResponse_2.corr, level);
    }
    psychoJS.experiment.addData('confirmResponse_2.keys', confirmResponse_2.keys);
    if (typeof confirmResponse_2.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('confirmResponse_2.rt', confirmResponse_2.rt);
        psychoJS.experiment.addData('confirmResponse_2.duration', confirmResponse_2.duration);
        }
    
    confirmResponse_2.stop();
    // the Routine "catchQuestion" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var _endingBreakKey_allKeys;
var runIntervalBreakComponents;
function runIntervalBreakRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'runIntervalBreak' ---
    t = 0;
    runIntervalBreakClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    routineTimer.add(60.000000);
    // update component parameters for each repeat
    psychoJS.experiment.addData('runIntervalBreak.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_3
    nAudioRun = (nAudioRun + 1);
    
    endingBreakKey.keys = undefined;
    endingBreakKey.rt = undefined;
    _endingBreakKey_allKeys = [];
    // keep track of which components have finished
    runIntervalBreakComponents = [];
    runIntervalBreakComponents.push(breakText_audio);
    runIntervalBreakComponents.push(endingBreakKey);
    
    for (const thisComponent of runIntervalBreakComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function runIntervalBreakRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'runIntervalBreak' ---
    // get current time
    t = runIntervalBreakClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *breakText_audio* updates
    if (t >= 0.0 && breakText_audio.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      breakText_audio.tStart = t;  // (not accounting for frame time here)
      breakText_audio.frameNStart = frameN;  // exact frame index
      
      breakText_audio.setAutoDraw(true);
    }
    
    frameRemains = 0.0 + 60 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (breakText_audio.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      breakText_audio.setAutoDraw(false);
    }
    
    // *endingBreakKey* updates
    if (t >= 0.5 && endingBreakKey.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      endingBreakKey.tStart = t;  // (not accounting for frame time here)
      endingBreakKey.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { endingBreakKey.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { endingBreakKey.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { endingBreakKey.clearEvents(); });
    }
    
    frameRemains = 60  - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if ((endingBreakKey.status === PsychoJS.Status.STARTED || endingBreakKey.status === PsychoJS.Status.FINISHED) && t >= frameRemains) {
      endingBreakKey.status = PsychoJS.Status.FINISHED;
        }
      
    if (endingBreakKey.status === PsychoJS.Status.STARTED) {
      let theseKeys = endingBreakKey.getKeys({keyList: ['space'], waitRelease: false});
      _endingBreakKey_allKeys = _endingBreakKey_allKeys.concat(theseKeys);
      if (_endingBreakKey_allKeys.length > 0) {
        endingBreakKey.keys = _endingBreakKey_allKeys[_endingBreakKey_allKeys.length - 1].name;  // just the last key pressed
        endingBreakKey.rt = _endingBreakKey_allKeys[_endingBreakKey_allKeys.length - 1].rt;
        endingBreakKey.duration = _endingBreakKey_allKeys[_endingBreakKey_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of runIntervalBreakComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine && routineTimer.getTime() > 0) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function runIntervalBreakRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'runIntervalBreak' ---
    for (const thisComponent of runIntervalBreakComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('runIntervalBreak.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(endingBreakKey.corr, level);
    }
    psychoJS.experiment.addData('endingBreakKey.keys', endingBreakKey.keys);
    if (typeof endingBreakKey.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('endingBreakKey.rt', endingBreakKey.rt);
        psychoJS.experiment.addData('endingBreakKey.duration', endingBreakKey.duration);
        routineTimer.reset();
        }
    
    endingBreakKey.stop();
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var _connectionEnding_allKeys;
var connectionComponents;
function connectionRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'connection' ---
    t = 0;
    connectionClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('connection.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_9
    connectionInstruct.setAlignHoriz('left');
    nAudioRun = 1;
    
    connectionEnding.keys = undefined;
    connectionEnding.rt = undefined;
    _connectionEnding_allKeys = [];
    // keep track of which components have finished
    connectionComponents = [];
    connectionComponents.push(connectionInstruct);
    connectionComponents.push(connectionEnding);
    
    for (const thisComponent of connectionComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function connectionRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'connection' ---
    // get current time
    t = connectionClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *connectionInstruct* updates
    if (t >= 0.0 && connectionInstruct.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      connectionInstruct.tStart = t;  // (not accounting for frame time here)
      connectionInstruct.frameNStart = frameN;  // exact frame index
      
      connectionInstruct.setAutoDraw(true);
    }
    
    
    // *connectionEnding* updates
    if (t >= 0.3 && connectionEnding.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      connectionEnding.tStart = t;  // (not accounting for frame time here)
      connectionEnding.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { connectionEnding.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { connectionEnding.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { connectionEnding.clearEvents(); });
    }
    
    if (connectionEnding.status === PsychoJS.Status.STARTED) {
      let theseKeys = connectionEnding.getKeys({keyList: ['space'], waitRelease: false});
      _connectionEnding_allKeys = _connectionEnding_allKeys.concat(theseKeys);
      if (_connectionEnding_allKeys.length > 0) {
        connectionEnding.keys = _connectionEnding_allKeys[_connectionEnding_allKeys.length - 1].name;  // just the last key pressed
        connectionEnding.rt = _connectionEnding_allKeys[_connectionEnding_allKeys.length - 1].rt;
        connectionEnding.duration = _connectionEnding_allKeys[_connectionEnding_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of connectionComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function connectionRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'connection' ---
    for (const thisComponent of connectionComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('connection.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(connectionEnding.corr, level);
    }
    psychoJS.experiment.addData('connectionEnding.keys', connectionEnding.keys);
    if (typeof connectionEnding.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('connectionEnding.rt', connectionEnding.rt);
        psychoJS.experiment.addData('connectionEnding.duration', connectionEnding.duration);
        routineTimer.reset();
        }
    
    connectionEnding.stop();
    // the Routine "connection" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var readyForTextComponents;
function readyForTextRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'readyForText' ---
    t = 0;
    readyForTextClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    routineTimer.add(2.580000);
    // update component parameters for each repeat
    psychoJS.experiment.addData('readyForText.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_7
    nTrial = 0;
    catchQ_presented = false;
    
    // keep track of which components have finished
    readyForTextComponents = [];
    readyForTextComponents.push(readyText_text);
    
    for (const thisComponent of readyForTextComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function readyForTextRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'readyForText' ---
    // get current time
    t = readyForTextClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *readyText_text* updates
    if (t >= 0.0 && readyText_text.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      readyText_text.tStart = t;  // (not accounting for frame time here)
      readyText_text.frameNStart = frameN;  // exact frame index
      
      readyText_text.setAutoDraw(true);
    }
    
    frameRemains = 0.0 + 2.58 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (readyText_text.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      readyText_text.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of readyForTextComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine && routineTimer.getTime() > 0) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function readyForTextRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'readyForText' ---
    for (const thisComponent of readyForTextComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('readyForText.stopped', globalClock.getTime());
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var textToPre;
var nWords;
var nScreens;
var currentScreen;
var tChange;
var startWord;
var textsOnScreen;
var tEnd;
var resetBool;
var textPreComponents;
var mouseInSlider;

function textPreRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'textPre' ---
    t = 0;
    textPreClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('textPre.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_8
    if ((((parID % 2) === 0) && (nTextRun === 1))) {
        textToPre = text1_1;
        psychoJS.experiment.addData('run_dup', 3);
    } else {
        if ((((parID % 2) === 1) && (nTextRun === 1))) {
            textToPre = text2_1;
            psychoJS.experiment.addData('run_dup', 3);
        } else {
            if ((((parID % 2) === 0) && (nTextRun === 2))) {
                textToPre = text1_2;
                psychoJS.experiment.addData('run_dup', 4);
            } else {
                if ((((parID % 2) === 1) && (nTextRun === 2))) {
                    textToPre = text2_2;
                    psychoJS.experiment.addData('run_dup', 4);
                }
            }
        }
    }
    textToPre = textToPre.split(" ");
    nWords = textToPre.length;
    nScreens = (Number.parseInt(((nWords - 1) / 10)) + 1);
    currentScreen = 0;
    tChange = (textPreClock.getTime() + 3.34);
    startWord = 0;
    textsOnScreen = textToPre.slice(startWord, (startWord + 10)).join(" ");
    tEnd = 9999999999999999999;
    resetBool = [];
    for (var i, _pj_c = 0, _pj_a = util.range((nScreens + 1)), _pj_b = _pj_a.length; (_pj_c < _pj_b); _pj_c += 1) {
        i = _pj_a[_pj_c];
        resetBool.push(0);
    }
    resetBool[currentScreen] = 1;
    currentScreen = (currentScreen + 1);
    thisFrame = 0;
    slider_data = [];
    continuousSlider_2.markerPos = 50;
    mouseRec = mouse.getPos();
    mouseInSlider = false;

    continuousSlider_2.reset()
    questionExplanation_2.setText(continuousQExp);
    // keep track of which components have finished
    textPreComponents = [];
    textPreComponents.push(texts);
    textPreComponents.push(slider_shape_2);
    textPreComponents.push(continuousSlider_2);
    textPreComponents.push(questionText_2);
    textPreComponents.push(questionExplanation_2);
    
    for (const thisComponent of textPreComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function textPreRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'textPre' ---
    // get current time
    t = textPreClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // Run 'Each Frame' code from code_8
    if (((nTextRun === 2) && (nTrial === 1))) {
      continueRoutine = false;
    }
    if (((! slider_shape_2.contains(mouse)) && (! mouseInSlider))) {
      if (((thisFrame % slideSpeed) === 0)) {
          slider_data.push([999, Number.parseInt((t * 1000))]);
      }
    }
    if ((slider_shape_2.contains(mouse) && (! mouseInSlider))) {
        mouseInSlider = true;
    }
    if ((slider_shape_2.contains(mouse) && mouseInSlider)) {
        if ((mouse.getPos()[slider_orientation] !== mouseRec[slider_orientation])) {
            mouseRec = mouse.getPos();
            continuousSlider_2.markerPos = (((mouseRec[slider_orientation] / slider_width) * 100) + (100 / 2));
        }
        if (continuousSlider_2.markerPos) {
            if ((oldRating !== continuousSlider_2.markerPos)) {
                oldRating = continuousSlider_2.markerPos;
            }
        }
        if ((continuousSlider_2.markerPos && ((thisFrame % slideSpeed) === 0))) {
            slider_data.push([util.round(oldRating, slider_decimals), Number.parseInt((t * 1000))]);
        }
    }
    thisFrame = (thisFrame + 1);
    if (((t >= (tChange - frameTolerance)) && (! resetBool[currentScreen]))) {
        resetBool[currentScreen] = 1;
        currentScreen = (currentScreen + 1);
        startWord = (startWord + 10);
        if (((startWord + 10) < nWords)) {
            tChange = (t + (10 / 3));
            textsOnScreen = textToPre.slice(startWord, (startWord + 10)).join(" ");
        } else {
            tChange = 99999999999999;
            tEnd = ((t + ((nWords - startWord) / 3)) + 1.5);
            textsOnScreen = textToPre.slice(startWord).join(" ");
        }
    }
    if ((t >= (tEnd - frameTolerance))) {
        continueRoutine = false;
    }
    
    
    if (texts.status === PsychoJS.Status.STARTED){ // only update if being drawn
      texts.setText(textsOnScreen, false);
    }
    
    // *texts* updates
    if (t >= 0.0 && texts.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      texts.tStart = t;  // (not accounting for frame time here)
      texts.frameNStart = frameN;  // exact frame index
      
      texts.setAutoDraw(true);
    }
    
    
    // *slider_shape_2* updates
    if (t >= 0.0 && slider_shape_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      slider_shape_2.tStart = t;  // (not accounting for frame time here)
      slider_shape_2.frameNStart = frameN;  // exact frame index
      
      slider_shape_2.setAutoDraw(true);
    }
    
    
    // *continuousSlider_2* updates
    if (t >= 0.0 && continuousSlider_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      continuousSlider_2.tStart = t;  // (not accounting for frame time here)
      continuousSlider_2.frameNStart = frameN;  // exact frame index
      
      continuousSlider_2.setAutoDraw(true);
    }
    
    
    // *questionText_2* updates
    if (t >= 0.0 && questionText_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      questionText_2.tStart = t;  // (not accounting for frame time here)
      questionText_2.frameNStart = frameN;  // exact frame index
      
      questionText_2.setAutoDraw(true);
    }
    
    
    // *questionExplanation_2* updates
    if (t >= 0.0 && questionExplanation_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      questionExplanation_2.tStart = t;  // (not accounting for frame time here)
      questionExplanation_2.frameNStart = frameN;  // exact frame index
      
      questionExplanation_2.setAutoDraw(true);
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of textPreComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function textPreRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'textPre' ---
    for (const thisComponent of textPreComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('textPre.stopped', globalClock.getTime());
    // Run 'End Routine' code from code_8
    psychoJS.experiment.addData("continuousResponses_text", slider_data);
    
    // the Routine "textPre" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var _endingBreakKey_2_allKeys;
var runBreakTextComponents;
function runBreakTextRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'runBreakText' ---
    t = 0;
    runBreakTextClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    routineTimer.add(60.000000);
    // update component parameters for each repeat
    psychoJS.experiment.addData('runBreakText.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_10
    nTextRun = (nTextRun + 1);
    nAudioRun = (nAudioRun + 1);
    
    endingBreakKey_2.keys = undefined;
    endingBreakKey_2.rt = undefined;
    _endingBreakKey_2_allKeys = [];
    // keep track of which components have finished
    runBreakTextComponents = [];
    runBreakTextComponents.push(breakText_text);
    runBreakTextComponents.push(endingBreakKey_2);
    
    for (const thisComponent of runBreakTextComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function runBreakTextRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'runBreakText' ---
    // get current time
    t = runBreakTextClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // Run 'Each Frame' code from code_10
    if ((nTextRun === 3)) {
        continueRoutine = false;
    }
    
    
    // *breakText_text* updates
    if (t >= 0.0 && breakText_text.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      breakText_text.tStart = t;  // (not accounting for frame time here)
      breakText_text.frameNStart = frameN;  // exact frame index
      
      breakText_text.setAutoDraw(true);
    }
    
    frameRemains = 0.0 + 60 - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (breakText_text.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      breakText_text.setAutoDraw(false);
    }
    
    // *endingBreakKey_2* updates
    if (t >= 0.5 && endingBreakKey_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      endingBreakKey_2.tStart = t;  // (not accounting for frame time here)
      endingBreakKey_2.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { endingBreakKey_2.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { endingBreakKey_2.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { endingBreakKey_2.clearEvents(); });
    }
    
    frameRemains = 60  - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if ((endingBreakKey_2.status === PsychoJS.Status.STARTED || endingBreakKey_2.status === PsychoJS.Status.FINISHED) && t >= frameRemains) {
      endingBreakKey_2.status = PsychoJS.Status.FINISHED;
        }
      
    if (endingBreakKey_2.status === PsychoJS.Status.STARTED) {
      let theseKeys = endingBreakKey_2.getKeys({keyList: ['space'], waitRelease: false});
      _endingBreakKey_2_allKeys = _endingBreakKey_2_allKeys.concat(theseKeys);
      if (_endingBreakKey_2_allKeys.length > 0) {
        endingBreakKey_2.keys = _endingBreakKey_2_allKeys[_endingBreakKey_2_allKeys.length - 1].name;  // just the last key pressed
        endingBreakKey_2.rt = _endingBreakKey_2_allKeys[_endingBreakKey_2_allKeys.length - 1].rt;
        endingBreakKey_2.duration = _endingBreakKey_2_allKeys[_endingBreakKey_2_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of runBreakTextComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine && routineTimer.getTime() > 0) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function runBreakTextRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'runBreakText' ---
    for (const thisComponent of runBreakTextComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('runBreakText.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(endingBreakKey_2.corr, level);
    }
    psychoJS.experiment.addData('endingBreakKey_2.keys', endingBreakKey_2.keys);
    if (typeof endingBreakKey_2.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('endingBreakKey_2.rt', endingBreakKey_2.rt);
        psychoJS.experiment.addData('endingBreakKey_2.duration', endingBreakKey_2.duration);
        routineTimer.reset();
        }
    
    endingBreakKey_2.stop();
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


var _endingEnding_allKeys;
var endingComponents;
function endingRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //--- Prepare to start Routine 'ending' ---
    t = 0;
    endingClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    psychoJS.experiment.addData('ending.started', globalClock.getTime());
    // Run 'Begin Routine' code from code_11
    endingText.setAlignHoriz('left');
    
    endingEnding.keys = undefined;
    endingEnding.rt = undefined;
    _endingEnding_allKeys = [];
    // keep track of which components have finished
    endingComponents = [];
    endingComponents.push(endingText);
    endingComponents.push(endingEnding);
    
    for (const thisComponent of endingComponents)
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
    return Scheduler.Event.NEXT;
  }
}


function endingRoutineEachFrame() {
  return async function () {
    //--- Loop for each frame of Routine 'ending' ---
    // get current time
    t = endingClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *endingText* updates
    if (t >= 0.0 && endingText.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      endingText.tStart = t;  // (not accounting for frame time here)
      endingText.frameNStart = frameN;  // exact frame index
      
      endingText.setAutoDraw(true);
    }
    
    
    // *endingEnding* updates
    if (t >= 0.3 && endingEnding.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      endingEnding.tStart = t;  // (not accounting for frame time here)
      endingEnding.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { endingEnding.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { endingEnding.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { endingEnding.clearEvents(); });
    }
    
    if (endingEnding.status === PsychoJS.Status.STARTED) {
      let theseKeys = endingEnding.getKeys({keyList: ['space'], waitRelease: false});
      _endingEnding_allKeys = _endingEnding_allKeys.concat(theseKeys);
      if (_endingEnding_allKeys.length > 0) {
        endingEnding.keys = _endingEnding_allKeys[_endingEnding_allKeys.length - 1].name;  // just the last key pressed
        endingEnding.rt = _endingEnding_allKeys[_endingEnding_allKeys.length - 1].rt;
        endingEnding.duration = _endingEnding_allKeys[_endingEnding_allKeys.length - 1].duration;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    for (const thisComponent of endingComponents)
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
        break;
      }
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function endingRoutineEnd(snapshot) {
  return async function () {
    //--- Ending Routine 'ending' ---
    for (const thisComponent of endingComponents) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    }
    psychoJS.experiment.addData('ending.stopped', globalClock.getTime());
    // update the trial handler
    if (currentLoop instanceof MultiStairHandler) {
      currentLoop.addResponse(endingEnding.corr, level);
    }
    psychoJS.experiment.addData('endingEnding.keys', endingEnding.keys);
    if (typeof endingEnding.keys !== 'undefined') {  // we had a response
        psychoJS.experiment.addData('endingEnding.rt', endingEnding.rt);
        psychoJS.experiment.addData('endingEnding.duration', endingEnding.duration);
        routineTimer.reset();
        }
    
    endingEnding.stop();
    // the Routine "ending" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    // Routines running outside a loop should always advance the datafile row
    if (currentLoop === psychoJS.experiment) {
      psychoJS.experiment.nextEntry(snapshot);
    }
    return Scheduler.Event.NEXT;
  }
}


function importConditions(currentLoop) {
  return async function () {
    psychoJS.importAttributes(currentLoop.getCurrentTrial());
    return Scheduler.Event.NEXT;
    };
}


async function quitPsychoJS(message, isCompleted) {
  // Check for and save orphaned data
  if (psychoJS.experiment.isEntryEmpty()) {
    psychoJS.experiment.nextEntry();
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  psychoJS.window.close();
  psychoJS.quit({message: message, isCompleted: isCompleted});
  
  return Scheduler.Event.QUIT;
}
